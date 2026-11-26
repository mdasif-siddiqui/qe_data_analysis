"""Job pipeline: SCF -> bands -> bands.x -> NSCF -> dos.x -> projwfc.x, plus convergence tests.

One job runs at a time (a queue), so calculations never compete for cores.
Every job lives in runs/<id>/ with its inputs, outputs, job.json (status) and results.json.
"""
from __future__ import annotations

import json
import math
import os
import queue
import re
import shutil
import signal
import subprocess
import threading
import time
from datetime import datetime
from pathlib import Path

import numpy as np

from . import cif, config, nmr, optics, phonon, pseudos, qein, qeout, spectra

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "runs"
# QE's scratch (outdir: wavefunctions, charge density) lives on the Linux filesystem: reading and
# writing it under /mnt/c is several times slower. Inputs/outputs stay in runs/<id>/.
SCRATCH = Path.home() / ".cache" / "qe_studio"
QE_BIN = Path(config.QE_BIN) if config.QE_BIN else Path.home() / "miniforge3" / "envs" / "qe" / "bin"
if not (QE_BIN / "pw.x").exists() and shutil.which("pw.x"):
    QE_BIN = Path(shutil.which("pw.x")).parent


class LimitError(ValueError):
    """Input exceeds the configured limits (public mode)."""


def _cpu_quota() -> int | None:
    """CPUs allowed by the container (cgroup v2 cpu.max), if limited."""
    try:
        q, period = open("/sys/fs/cgroup/cpu.max").read().split()
        if q != "max":
            return max(1, int(int(q) / int(period)))
    except (OSError, ValueError):
        pass
    return None


def physical_cores() -> int:
    n = None
    try:
        out = subprocess.run(["lscpu", "-p=CORE,SOCKET"], capture_output=True, text=True).stdout
        cores = {l for l in out.splitlines() if l and not l.startswith("#")}
        if cores:
            n = len(cores)
    except OSError:
        pass
    if n is None:
        n = max(1, (os.cpu_count() or 2) // 2)
    try:
        n = min(n, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        pass
    q = _cpu_quota()
    if q:
        n = min(n, q)
    if config.MAX_NP:
        n = min(n, config.MAX_NP)
    return max(1, n)


def exe(name: str) -> str:
    p = QE_BIN / name
    return str(p) if p.exists() else (shutil.which(name) or name)


def qe_version() -> str | None:
    try:
        import tempfile
        with tempfile.TemporaryDirectory() as td:  # pw.x with empty input leaves CRASH/input_tmp.in behind
            r = subprocess.run([exe("pw.x")], input="", capture_output=True, text=True, timeout=60,
                               env={**os.environ, "OMP_NUM_THREADS": "1"}, cwd=td)
        m = re.search(r"Program PWSCF v\.?\s*([\w.]+)", r.stdout)
        return m.group(1) if m else None
    except Exception:
        return None


# --------------------------------------------------------------------------- job object
class Job:
    def __init__(self, jid: str):
        self.id = jid
        self.dir = RUNS / jid
        self.scratch = SCRATCH / jid
        self.proc: subprocess.Popen | None = None
        self.cancelled = False

    # persistent state -------------------------------------------------------
    @property
    def state_file(self):
        return self.dir / "job.json"

    def load(self) -> dict:
        for _ in range(5):
            try:
                return json.loads(self.state_file.read_text())
            except (OSError, json.JSONDecodeError):
                time.sleep(0.05)
        return {}

    def save(self, st: dict):
        tmp = self.state_file.with_suffix(".tmp")
        tmp.write_text(json.dumps(st, indent=1))
        tmp.replace(self.state_file)

    def update(self, **kw):
        st = self.load()
        st.update(kw)
        self.save(st)
        return st

    def step(self, name: str, **kw):
        st = self.load()
        for s in st.get("steps", []):
            if s["name"] == name:
                s.update(kw)
        self.save(st)

    def log(self, msg: str):
        with open(self.dir / "pipeline.log", "a") as f:
            f.write(f"[{datetime.now():%H:%M:%S}] {msg}\n")


# --------------------------------------------------------------------------- manager
class Manager:
    def __init__(self):
        RUNS.mkdir(exist_ok=True)
        self.q: "queue.Queue[str]" = queue.Queue()
        self.jobs: dict[str, Job] = {}
        self.current: Job | None = None
        self.launcher: str | None = None  # remembered working launch mode
        # jobs interrupted by a restart are marked as such
        for d in RUNS.iterdir():
            j = Job(d.name)
            st = j.load()
            if st.get("status") in ("running", "queued"):
                j.update(status="failed", error="Interrupted: the server was restarted.")
        threading.Thread(target=self._worker, daemon=True).start()

    def get(self, jid: str) -> Job:
        if jid not in self.jobs:
            self.jobs[jid] = Job(jid)
        return self.jobs[jid]

    def submit(self, job: Job):
        self.jobs[job.id] = job
        self.q.put(job.id)

    def cancel(self, jid: str):
        job = self.get(jid)
        job.cancelled = True
        st = job.load()
        if st.get("status") == "queued":
            job.update(status="cancelled")
        if job.proc and job.proc.poll() is None:
            try:
                os.killpg(job.proc.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass

    def _worker(self):
        while True:
            jid = self.q.get()
            job = self.get(jid)
            if job.cancelled or job.load().get("status") == "cancelled":
                continue
            self.current = job
            try:
                run_pipeline(self, job)
            except Exception as e:  # report, never kill the worker
                job.log(f"ERROR: {e}")
                job.update(status="cancelled" if job.cancelled else "failed", error=str(e),
                           finished=time.time())
            finally:
                self.current = None

    # process execution ---------------------------------------------------------
    def run(self, job: Job, program: str, inp: Path, out: Path, np_: int, use_mpi: bool = True,
            cwd: Path | None = None, gipaw: bool = False) -> int:
        """Run a QE program, input via -in (pw.x, gipaw.x) or stdin (post-processing).

        gipaw=True runs the self-built pw.x/gipaw.x from ~/gipaw/bin with their own MPI (build_gipaw.sh)."""
        cwd = cwd or job.dir
        env = {**os.environ, "OMP_NUM_THREADS": "1",
               "PATH": f"{QE_BIN}:{os.environ.get('PATH', '')}"}
        if gipaw:
            env.update(nmr.env())
        # launch modes, tried in order until the program actually starts:
        #   "mpi" -> mpirun -np N ; "mpi_cma" -> same with CMA shared memory disabled ; "serial"
        if use_mpi and np_ > 1:
            modes = [self.launcher] if self.launcher else ["mpi", "mpi_cma", "serial"]
        else:
            modes = ["serial"]
        rc = 1
        for mode in modes:
            env_m = dict(env)
            launcher: list[str] = []
            if mode in ("mpi", "mpi_cma"):
                launcher = [str(nmr.BUILD_ENV / "bin" / "mpirun") if gipaw else exe("mpirun"), "-np", str(np_)]
                if os.geteuid() == 0:
                    launcher.insert(1, "--allow-run-as-root")
            if mode == "mpi_cma":
                env_m["OMPI_MCA_smsc"] = "^cma"
                env_m["OMPI_MCA_btl_vader_single_copy_mechanism"] = "none"
            cmd = launcher + [str(nmr.GIPAW_BIN / program) if gipaw else exe(program)]
            stdin = subprocess.DEVNULL
            if program in ("pw.x", "gipaw.x"):
                cmd += ["-in", inp.name]
            else:
                stdin = open(inp)
            job.log(f"$ {'OMP_NUM_THREADS=1 ' + ' '.join(cmd)}{'' if program in ('pw.x', 'gipaw.x') else ' < ' + inp.name} > {out.name}")
            with open(out, "w") as fo, open(out.with_suffix(".err"), "w") as fe:
                job.proc = subprocess.Popen(cmd, cwd=cwd, stdout=fo, stderr=fe, stdin=stdin, env=env_m,
                                            start_new_session=True)
                try:
                    rc = job.proc.wait(timeout=config.MAX_STEP_SECONDS or None)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(job.proc.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    job.proc.wait()
                    raise RuntimeError(f"{program} exceeded the time limit of {config.MAX_STEP_SECONDS} s for one step.")
            if stdin is not subprocess.DEVNULL:
                stdin.close()
            if job.cancelled:
                raise RuntimeError("Cancelled by user")
            text = out.read_text(errors="replace")
            started = re.search(r"Program \S+ v\.", text) is not None
            if started or mode == "serial":
                if use_mpi and np_ > 1 and self.launcher is None and started:
                    self.launcher = mode
                    job.log(f"Launch mode '{mode}' works; using it for later runs.")
                return rc
            err = out.with_suffix(".err").read_text(errors="replace")[-600:]
            job.log(f"MPI launch failed ({err.strip().splitlines()[-1] if err.strip() else 'no output'}); trying next mode.")
        return rc


MANAGER: Manager | None = None


def manager() -> Manager:
    global MANAGER
    if MANAGER is None:
        MANAGER = Manager()
    return MANAGER


# --------------------------------------------------------------------------- creation
def create_job(scf_text: str, scf_name: str, pseudo_files: dict[str, bytes], options: dict) -> Job:
    mgr = manager()  # start the manager before writing a 'queued' job (it resets stale ones on start)
    fixed, notes = qein.normalize_text(scf_text)
    parsed = qein.parse_input(fixed)
    check_limits(fixed, parsed, options)
    for sp in qein.species_list(parsed):
        pseudos.safe_name(sp["pseudo"])  # raises ValueError for names with paths
    if mgr.q.qsize() >= config.MAX_QUEUE:
        raise LimitError(f"The queue is full ({config.MAX_QUEUE} jobs waiting). Try again later.")
    cleanup_old_runs()
    prefix = str(qein.get(parsed, "control", "prefix", "pwscf"))
    safe = re.sub(r"[^\w.-]", "_", prefix)[:40]
    base = f"{datetime.now():%Y%m%d-%H%M%S}_{safe}"
    if config.PUBLIC:  # unguessable ids: others cannot open your run by guessing its time stamp
        import secrets
        base += "_" + secrets.token_hex(6)
    jid, n = base, 1
    while (RUNS / jid).exists():
        n += 1
        jid = f"{base}-{n}"
    job = Job(jid)
    job.dir.mkdir(parents=True)
    (job.dir / ("original_" + os.path.basename(scf_name or "scf.in"))).write_text(scf_text)
    pdir = job.dir / "pseudos"
    pdir.mkdir()
    for name, data in pseudo_files.items():
        (pdir / pseudos.safe_name(name)).write_bytes(data)
    steps = [{"name": "prepare", "label": "Prepare inputs & pseudopotentials"}]
    rmode = relax_mode(fixed, options)
    if rmode:
        steps += [{"name": "relax", "label": f"Structure relaxation (pw.x {rmode})"}]
    steps += [{"name": "scf", "label": "SCF (pw.x)" + (" on the relaxed structure" if rmode else "")}]
    if options.get("bands", True):
        steps += [{"name": "bands", "label": "Band structure (pw.x bands)"},
                  {"name": "bandsx", "label": "bands.x post-processing"}]
    if options.get("dos", True):
        steps += [{"name": "nscf", "label": "NSCF on dense k-grid (pw.x)"},
                  {"name": "dos", "label": "Density of states (dos.x)"},
                  {"name": "pdos", "label": "Projected DOS (projwfc.x)"}]
    if options.get("phonons"):
        steps += [{"name": "ph_scf", "label": "Phonons: tight SCF (pw.x)"},
                  {"name": "ph", "label": "Phonons: DFPT on q-grid (ph.x)"},
                  {"name": "q2r", "label": "Force constants (q2r.x)"},
                  {"name": "matdyn", "label": "Dispersion + phonon DOS (matdyn.x)"}]
    if options.get("spectra"):
        steps += [{"name": "sp_scf", "label": "IR/Raman: tight SCF (pw.x)"},
                  {"name": "sp_ph", "label": "IR/Raman: Γ phonons, dielectric, Raman (ph.x)"},
                  {"name": "dynmat", "label": "IR/Raman intensities (dynmat.x)"}]
    if options.get("uvvis"):
        steps += [{"name": "uv_nscf", "label": "UV-Vis: SCF + NSCF with empty bands (pw.x)"},
                  {"name": "optics", "label": "UV-Vis: dielectric function and absorption"}]
    if options.get("nmr"):
        if not nmr.available() and not config.PUBLIC:
            steps += [{"name": "gipaw_build", "label": "NMR: one-time GIPAW build (30-60 min, first NMR run only)"}]
        steps += [{"name": "nmr_scf", "label": "NMR: SCF (pw.x, GIPAW build)"},
                  {"name": "nmr", "label": "NMR shielding tensors (gipaw.x)"}]
    if options.get("conv"):
        steps += [{"name": "conv_ecut", "label": "Cutoff convergence test"},
                  {"name": "conv_k", "label": "k-point convergence test"}]
    steps.append({"name": "analysis", "label": "Analysis"})
    for s in steps:
        s["status"] = "pending"
    job.save({"id": jid, "prefix": prefix, "input_name": scf_name, "status": "queued",
              "created": time.time(), "options": options, "steps": steps, "notes": notes, "relax_mode": rmode})
    (job.dir / "scf.in").write_text(fixed)
    job.log(f"Job created from {scf_name}. " + " ".join(notes))
    mgr.submit(job)
    return job


def check_limits(text: str, parsed: dict, options: dict):
    """Raise LimitError if the job is bigger than this server allows (see config.py)."""
    sysnl = parsed["namelists"].get("system", {})
    nat = int(sysnl.get("nat", 0) or 0)
    if nat > config.MAX_ATOMS:
        raise LimitError(f"nat = {nat}: this server allows at most {config.MAX_ATOMS} atoms.")
    ecut = sysnl.get("ecutwfc")
    if isinstance(ecut, (int, float)) and ecut > config.MAX_ECUTWFC:
        raise LimitError(f"ecutwfc = {ecut} Ry: this server allows at most {config.MAX_ECUTWFC} Ry.")
    grids = []
    k = qein.kpoints_info(parsed)
    if k.get("grid"):
        grids.append(k["grid"])
    if options.get("nscf_grid"):
        grids.append([int(x) for x in options["nscf_grid"]])
    for g in grids:
        if g and int(np.prod(g)) > config.MAX_KPOINTS:
            raise LimitError(f"k-grid {'x'.join(map(str, g))} is larger than this server allows ({config.MAX_KPOINTS} points).")
    if options.get("phonons") and options.get("ph_qgrid"):
        qg = [int(x) for x in options["ph_qgrid"]]
        if config.PUBLIC and int(np.prod(qg)) > 64:
            raise LimitError("Phonon q-grid too large for this server (at most 4x4x4).")
    if config.PUBLIC:
        for key in ("ecut_list", "k_list"):
            vals = [v for v in re.split(r"[,\s]+", str(options.get(key) or "")) if v]
            if len(vals) > 8:
                raise LimitError("At most 8 values per convergence test.")
            if key == "ecut_list" and any(float(v) > config.MAX_ECUTWFC for v in vals):
                raise LimitError(f"Convergence cutoffs must be at most {config.MAX_ECUTWFC} Ry.")
            if key == "k_list" and any(int(float(v)) ** 3 > config.MAX_KPOINTS for v in vals):
                raise LimitError("A convergence k-point value is too large for this server.")


def cleanup_old_runs():
    if not config.KEEP_DAYS or not RUNS.exists():
        return
    cutoff = time.time() - config.KEEP_DAYS * 86400
    for d in RUNS.iterdir():
        st = Job(d.name).load()
        if st.get("status") in ("running", "queued"):
            continue
        if (st.get("created") or d.stat().st_mtime) < cutoff:
            shutil.rmtree(d, ignore_errors=True)
            shutil.rmtree(SCRATCH / d.name, ignore_errors=True)


# --------------------------------------------------------------------------- relaxation
RELAX_MODES = ("relax", "vc-relax")


def relax_mode(text: str, opts: dict) -> str | None:
    """'relax' / 'vc-relax' when the user asked for it or the input itself is a relaxation, else None."""
    mode = str(opts.get("relax") or "").lower()
    if mode in RELAX_MODES:
        return mode
    calc = str(qein.get(qein.parse_input(text), "control", "calculation", "scf")).lower()
    return calc if calc in RELAX_MODES else None


def relax_input(text: str, mode: str, outdir: str) -> str:
    """pw.x relax/vc-relax input with thresholds tight enough for phonons (forces 1e-4 Ry/bohr, 0.5 kbar)."""
    parsed = qein.parse_input(text)
    t = qein.set_param(text, "control", "calculation", mode)
    t = qein.set_param(t, "control", "outdir", outdir)
    t = qein.set_param(t, "control", "tprnfor", True)
    t = qein.set_param(t, "control", "tstress", True)
    t = qein.set_param(t, "control", "forc_conv_thr", min(float(qein.get(parsed, "control", "forc_conv_thr", 1e-4)), 1e-4))
    t = qein.set_param(t, "control", "etot_conv_thr", min(float(qein.get(parsed, "control", "etot_conv_thr", 1e-5)), 1e-5))
    t = qein.set_param(t, "control", "nstep", max(int(qein.get(parsed, "control", "nstep", 100) or 100), 100))
    t = qein.set_param(t, "electrons", "conv_thr", min(float(qein.get(parsed, "electrons", "conv_thr", 1e-10)), 1e-10))
    if qein._namelist_span(t, "ions") is None:
        t = qein.set_param(t, "ions", "ion_dynamics", "bfgs")
    if mode == "vc-relax":
        t = qein.set_param(t, "cell", "press_conv_thr", min(float(qein.get(parsed, "cell", "press_conv_thr", 0.5)), 0.5))
    return t


def relaxed_input(text: str, out: Path, mode: str) -> tuple[str, dict]:
    """The scf input with the final structure of a relax/vc-relax run (last structure if it did not converge).

    Positions are written in crystal coordinates; for vc-relax the cell becomes ibrav = 0 with CELL_PARAMETERS in
    angstrom, so the result does not depend on how the original input defined alat. if_pos flags are kept."""
    from ase.io import read as ase_read
    out_text = out.read_text(errors="replace")
    final = ase_read(str(out), format="espresso-out", index=-1)
    before = qein.read_structure(text)
    parsed = qein.parse_input(text)
    t = qein.set_param(text, "control", "calculation", "scf")
    for nl, keys in (("control", ("forc_conv_thr", "etot_conv_thr", "nstep")), ("ions", ()), ("cell", ())):
        for k in keys:
            t = qein.remove_param(t, nl, k)
    cell = final.cell.array if mode == "vc-relax" else before.cell.array
    if mode == "vc-relax":
        t = qein.set_param(t, "system", "ibrav", 0)
        for k in ["A", "B", "C", "cosAB", "cosAC", "cosBC"] + [f"celldm({i})" for i in range(1, 7)]:
            t = qein.remove_param(t, "system", k)
        t = qein.replace_card(t, "CELL_PARAMETERS", "CELL_PARAMETERS angstrom\n" +
                              "\n".join("  " + " ".join(f"{x:15.10f}" for x in v) for v in cell) + "\n")
    frac = np.linalg.solve(cell.T, final.positions.T).T % 1.0
    frac[np.isclose(frac, 1.0, atol=1e-9)] = 0.0
    card = parsed["cards"]["ATOMIC_POSITIONS"]["lines"]
    lines = []
    for i, (lab, f) in enumerate(zip(before.info["labels"], frac)):
        flags = card[i].split()[4:7] if i < len(card) else []
        lines.append(f"  {lab:4s} " + " ".join(f"{x:14.10f}" for x in f) + ("  " + " ".join(flags) if len(flags) == 3 else ""))
    t = qein.replace_card(t, "ATOMIC_POSITIONS", "ATOMIC_POSITIONS crystal\n" + "\n".join(lines) + "\n")
    # the &IONS / &CELL namelists are not allowed for calculation = 'scf' in some QE versions: drop them
    for nl in ("ions", "cell"):
        sp = qein._namelist_span(t, nl)
        if sp:
            endline = t.find("\n", sp[2])
            t = t[:sp[0]] + t[(len(t) if endline < 0 else endline + 1):]
    p = re.findall(r"P=\s*(" + r"[-\d.]+" + r")", out_text)
    f = re.findall(r"Total force =\s*([-\d.]+)", out_text)
    summary = {"mode": mode, "converged": "End final coordinates" in out_text or "bfgs converged" in out_text,
               "bfgs_steps": len(re.findall(r"^!", out_text, re.M)),
               "volume_before_A3": float(before.get_volume()), "volume_after_A3": float(abs(np.linalg.det(cell))),
               "abc_before": [float(x) for x in before.cell.lengths()], "abc_after": [float(x) for x in np.linalg.norm(cell, axis=1)],
               "pressure_kbar_first": float(p[0]) if p else None, "pressure_kbar": float(p[-1]) if p else None,
               "force_first": float(f[0]) if f else None, "force": float(f[-1]) if f else None,
               "energy_change_ry": None}
    e = re.findall(r"^!\s+total energy\s+=\s+([-\d.]+)", out_text, re.M)
    if len(e) >= 2:
        summary["energy_change_ry"] = float(e[-1]) - float(e[0])
    return t, summary


# --------------------------------------------------------------------------- pipeline
def _check(job: Job, step: str, out: Path, rc: int, required=True) -> bool:
    text = out.read_text(errors="replace") if out.exists() else ""
    ok = "JOB DONE" in text
    errs = qeout.errors_of(text)
    if "convergence NOT achieved" in text:
        ok = False
        errs.append("SCF convergence NOT achieved (try smaller mixing_beta, more electron_maxstep, or smearing).")
    if ok:
        return True
    tail = "\n".join(text.splitlines()[-15:])
    msg = "; ".join(errs) if errs else f"exit code {rc}, no 'JOB DONE'. Last lines:\n{tail}"
    hint = _hint(text + (out.with_suffix(".err").read_text(errors="replace") if out.with_suffix(".err").exists() else ""))
    job.step(step, status="failed", error=msg + (f"\nHint: {hint}" if hint else ""))
    job.log(f"{step} FAILED: {msg}")
    if required:
        raise RuntimeError(f"{step} failed: {msg}" + (f" — {hint}" if hint else ""))
    return False


def _hint(text: str) -> str:
    rules = [
        ("read_rhog", "tmp/ was written by a different QE version or a failed run: rerun the whole job."),
        ("reading bands namelist", "bands.x input is malformed (missing newline after '/')."),
        ("not found", "a file is missing (pseudopotential or save folder)."),
        ("buffer overflow", "this is the broken apt QE 6.7 – make sure the conda QE is used."),
        ("not enough slots", "lower the number of MPI processes."),
        ("cholesky", "diagonalization problem: try diagonalization='cg' or a different mixing_beta."),
        ("S matrix not positive definite", "in pw.x: atoms may be too close (check positions / cell); in gipaw.x: "
                                           "try diagonalization = 'cg' (QE Studio retries with it automatically)."),
        ("too many bands are not converged", "try diagonalization='cg' or a larger ecutwfc."),
    ]
    for key, h in rules:
        if key.lower() in text.lower():
            return h
    return ""


def _timed(job: Job, name: str, fn):
    t0 = time.time()
    job.step(name, status="running", started=t0)
    job.update(current_step=name)
    job.log(f"--- {name} ---")
    res = fn()
    cur = next((s for s in job.load().get("steps", []) if s["name"] == name), {})
    if cur.get("status") in ("failed", "skipped"):
        job.step(name, seconds=round(time.time() - t0, 2))
    else:
        job.step(name, status="done", seconds=round(time.time() - t0, 2))
    return res


def run_pipeline(mgr: Manager, job: Job):
    st = job.update(status="running", started=time.time())
    opts = st["options"]
    np_ = max(1, int(opts.get("np") or physical_cores()))
    if config.PUBLIC:
        np_ = min(np_, physical_cores())
    d = job.dir
    scf_in = d / "scf.in"
    text = scf_in.read_text()
    parsed = qein.parse_input(text)
    prefix = str(qein.get(parsed, "control", "prefix", "pwscf"))
    results: dict = {"prefix": prefix, "np": np_}
    outdir = str(job.scratch) + "/"

    # ---- prepare
    def prepare():
        t = text
        t = qein.set_param(t, "control", "outdir", outdir)
        t = qein.set_param(t, "control", "pseudo_dir", "./pseudos/")
        t = qein.remove_param(t, "control", "wfcdir")  # would write outside the run folder
        if config.MAX_STEP_SECONDS:
            t = qein.set_param(t, "control", "max_seconds", float(max(60, config.MAX_STEP_SECONDS - 30)))
        calc = str(qein.get(parsed, "control", "calculation", "scf")).lower()
        if calc != "scf":       # relax/vc-relax run in their own step; scf.in then gets the relaxed structure
            t = qein.set_param(t, "control", "calculation", "scf")
        scf_in.write_text(t)
        names = [s["pseudo"] for s in qein.species_list(parsed)]
        extra = []
        orig_pd = qein.get(parsed, "control", "pseudo_dir")
        if orig_pd and os.path.isabs(str(orig_pd)):
            extra.append(Path(str(orig_pd)))
        info = pseudos.ensure(names, d / "pseudos", extra=extra,
                              allow_download=opts.get("download_pseudos", True), log=job.log)
        results["pseudos"] = info
        shutil.rmtree(job.scratch, ignore_errors=True)  # never reuse a save dir from another run/version
        job.scratch.mkdir(parents=True, exist_ok=True)
        job.log(f"Scratch (outdir) folder: {outdir}")
        job.log(f"Pseudopotentials ready: {', '.join(names)}")
    _timed(job, "prepare", prepare)

    # ---- relaxation: every later step (SCF, bands, DOS, phonons, IR/Raman, UV-Vis, NMR) uses the relaxed structure
    rmode = st.get("relax_mode")
    if rmode and any(s["name"] == "relax" for s in st["steps"]):
        def relax():
            rd = job.scratch / "relax"
            rd.mkdir(parents=True, exist_ok=True)
            (d / "relax.in").write_text(relax_input(scf_in.read_text(), rmode, str(rd) + "/"))
            rc = mgr.run(job, "pw.x", d / "relax.in", d / "relax.out", np_)
            out_text = (d / "relax.out").read_text(errors="replace") if (d / "relax.out").exists() else ""
            if not re.search(r"^!", out_text, re.M):
                _check(job, "relax", d / "relax.out", rc)          # nothing usable: report the error and stop
            new_text, summary = relaxed_input(scf_in.read_text(), d / "relax.out", rmode)
            scf_in.write_text(new_text)
            (d / "relax.json").write_text(json.dumps(summary, indent=1))
            shutil.rmtree(rd, ignore_errors=True)
            job.log(f"Relaxed structure ({rmode}, {summary['bfgs_steps']} SCF cycles, "
                    f"{'converged' if summary['converged'] else 'NOT converged - using the last structure'}): "
                    f"V {summary['volume_before_A3']:.3f} -> {summary['volume_after_A3']:.3f} A^3; scf.in updated.")
            if not summary["converged"]:
                job.step("relax", error="The relaxation did not converge (see relax.out); the last structure is used.")
        _timed(job, "relax", relax)

    # ---- SCF
    def scf():
        rc = mgr.run(job, "pw.x", scf_in, d / "scf.out", np_)
        _check(job, "scf", d / "scf.out", rc)
        if not (job.scratch / f"{prefix}.save" / "data-file-schema.xml").exists():
            raise RuntimeError(f"SCF finished but {outdir}{prefix}.save is missing.")
        return qeout.parse_pw((d / "scf.out").read_text())
    scf_res = _timed(job, "scf", scf)
    results["scf"] = _slim(scf_res)
    nelec = scf_res.get("nelec") or 0
    scf_text_final = scf_in.read_text()
    p_final = qein.parse_input(scf_text_final)
    occ = str(qein.get(p_final, "system", "occupations", "fixed")).lower()
    noncolin = bool(qein.get(p_final, "system", "noncolin", False))
    nocc = nelec if noncolin else nelec / 2
    nbnd_scf = scf_res.get("nbnd") or int(math.ceil(nocc))
    nbnd_extra = max(nbnd_scf, int(math.ceil(nocc)) + max(4, int(math.ceil(0.3 * nocc))))
    # For relax/vc-relax the later steps use the relaxed structure (QE reads it from the save dir).
    base_calc_text = scf_text_final

    # ---- bands
    if any(s["name"] == "bands" for s in st["steps"]):
        def bands():
            atoms = qein.read_structure(base_calc_text)
            path = (opts.get("kpath") or "").strip() or None
            kp = qein.make_kpath(atoms, path=path, density=float(opts.get("kdensity") or 40))
            t = qein.set_param(base_calc_text, "control", "calculation", "bands")
            t = qein.set_param(t, "control", "verbosity", "high")
            t = qein.set_param(t, "system", "nbnd", int(opts.get("nbnd_bands") or nbnd_extra))
            t = qein.replace_card(t, "K_POINTS", kp["card"])
            (d / "bands.in").write_text(t)
            (d / "kpath.json").write_text(json.dumps(kp, indent=1))
            rc = mgr.run(job, "pw.x", d / "bands.in", d / "bands.out", np_)
            _check(job, "bands", d / "bands.out", rc)
            return kp
        kp = _timed(job, "bands", bands)

        def bandsx():
            (d / "bands_pp.in").write_text(
                f"&BANDS\n  prefix = '{prefix}'\n  outdir = '{outdir}'\n  filband = '{prefix}_bands.dat'\n"
                f"  lsym = .false.\n/\n")
            rc = mgr.run(job, "bands.x", d / "bands_pp.in", d / "bands_pp.out", 1, use_mpi=False)
            _check(job, "bandsx", d / "bands_pp.out", rc, required=False)
        try:
            _timed(job, "bandsx", bandsx)
        except Exception as e:
            job.step("bandsx", status="failed", error=str(e))

    # ---- NSCF + DOS + PDOS
    if any(s["name"] == "nscf" for s in st["steps"]):
        def nscf():
            grid = opts.get("nscf_grid")
            if not grid:
                atoms = qein.read_structure(base_calc_text)
                kinfo = qein.kpoints_info(p_final)
                grid = [min(2 * g, 24) for g in kinfo["grid"]] if kinfo.get("grid") else qein.kgrid_for_density(atoms, 12)
            grid = [int(x) for x in grid]
            t = qein.set_param(base_calc_text, "control", "calculation", "nscf")
            t = qein.set_param(t, "control", "verbosity", "high")
            t = qein.set_param(t, "system", "nbnd", int(opts.get("nbnd_bands") or nbnd_extra))
            if occ != "smearing":
                t = qein.set_param(t, "system", "occupations", "tetrahedra_opt")
            t = qein.replace_card(t, "K_POINTS", "K_POINTS automatic\n  {} {} {} 0 0 0\n".format(*grid))
            (d / "nscf.in").write_text(t)
            rc = mgr.run(job, "pw.x", d / "nscf.in", d / "nscf.out", np_)
            _check(job, "nscf", d / "nscf.out", rc)
            results["nscf_grid"] = grid
        ok_nscf = True
        try:
            _timed(job, "nscf", nscf)
        except Exception as e:
            ok_nscf = False
            job.step("nscf", status="failed", error=str(e))
            for s in ("dos", "pdos"):
                job.step(s, status="skipped")
        if ok_nscf:
            def dos():
                (d / "dos.in").write_text(
                    f"&DOS\n  prefix = '{prefix}'\n  outdir = '{outdir}'\n  fildos = '{prefix}.dos'\n  DeltaE = 0.01\n/\n")
                rc = mgr.run(job, "dos.x", d / "dos.in", d / "dos.out", np_)
                _check(job, "dos", d / "dos.out", rc, required=False)
            def pdos():
                (d / "projwfc.in").write_text(
                    f"&PROJWFC\n  prefix = '{prefix}'\n  outdir = '{outdir}'\n  filpdos = '{prefix}'\n  DeltaE = 0.01\n/\n")
                rc = mgr.run(job, "projwfc.x", d / "projwfc.in", d / "projwfc.out", np_)
                _check(job, "pdos", d / "projwfc.out", rc, required=False)
            for name, fn in (("dos", dos), ("pdos", pdos)):
                try:
                    _timed(job, name, fn)
                except Exception as e:
                    job.step(name, status="failed", error=str(e))

    # ---- phonons
    if opts.get("phonons"):
        try:
            run_phonons(mgr, job, base_calc_text, prefix, opts, np_)
        except Exception as e:
            job.log(f"Phonons stopped: {e}")
            for s_ in ("ph_scf", "ph", "q2r", "matdyn"):
                cur = next((x for x in job.load()["steps"] if x["name"] == s_), {})
                if cur.get("status") in ("pending", "running"):
                    job.step(s_, status="skipped" if cur.get("status") == "pending" else "failed",
                             **({} if cur.get("status") == "pending" else {"error": str(e)}))

    # ---- IR / Raman
    if opts.get("spectra"):
        try:
            run_spectra(mgr, job, base_calc_text, prefix, scf_res, np_)
        except Exception as e:
            job.log(f"IR/Raman stopped: {e}")
            for s_ in ("sp_scf", "sp_ph", "dynmat"):
                cur = next((x for x in job.load()["steps"] if x["name"] == s_), {})
                if cur.get("status") in ("pending", "running"):
                    job.step(s_, status="skipped" if cur.get("status") == "pending" else "failed",
                             **({} if cur.get("status") == "pending" else {"error": str(e)}))

    # ---- UV-Vis
    if opts.get("uvvis"):
        try:
            run_uvvis(mgr, job, base_calc_text, prefix, opts, np_)
        except Exception as e:
            job.log(f"UV-Vis stopped: {e}")
            for s_ in ("uv_nscf", "optics"):
                cur = next((x for x in job.load()["steps"] if x["name"] == s_), {})
                if cur.get("status") in ("pending", "running"):
                    job.step(s_, status="skipped" if cur.get("status") == "pending" else "failed",
                             **({} if cur.get("status") == "pending" else {"error": str(e)}))

    # ---- NMR
    if opts.get("nmr"):
        try:
            run_nmr(mgr, job, base_calc_text, prefix, np_)
        except Exception as e:
            if job.cancelled:
                raise
            job.log(f"NMR stopped: {e}")
            for s_ in ("gipaw_build", "nmr_scf", "nmr"):
                cur = next((x for x in job.load()["steps"] if x["name"] == s_), {})
                if cur.get("status") in ("pending", "running"):
                    job.step(s_, status="skipped" if cur.get("status") == "pending" else "failed",
                             **({} if cur.get("status") == "pending" else {"error": str(e)}))

    # ---- convergence tests
    if opts.get("conv"):
        _timed(job, "conv_ecut", lambda: convergence(mgr, job, base_calc_text, "ecut", opts, np_))
        _timed(job, "conv_k", lambda: convergence(mgr, job, base_calc_text, "k", opts, np_))

    # ---- analysis
    _timed(job, "analysis", lambda: analyse_job(job))
    final = job.load()
    failed = [s["name"] for s in final["steps"] if s.get("status") == "failed"]
    job.update(status="done" if not failed else "done_with_errors", finished=time.time(), current_step=None,
               failed_steps=failed)
    if not opts.get("keep_scratch"):
        shutil.rmtree(job.scratch, ignore_errors=True)
    job.log("Pipeline finished" + (f" with failed steps: {failed}" if failed else "."))


def _slim(d: dict) -> dict:
    return {k: v for k, v in d.items() if k not in ("bands",)}


# --------------------------------------------------------------------------- phonons
def default_qgrid(atoms) -> list[int]:
    """About 4 q-points along the shortest reciprocal-space direction, fewer for big cells."""
    return qein.kgrid_for_density(atoms, 4 if len(atoms) <= 8 else 2)


def run_phonons(mgr: Manager, job: Job, text: str, prefix: str, opts: dict, np_: int):
    """Tight SCF -> ph.x (DFPT on a q-grid) -> q2r.x -> matdyn.x (path + DOS), in runs/<id>/phonon/.

    Uses its own scratch folder: ph.x needs the save directory of a plain SCF, and the bands/NSCF
    steps overwrite the wavefunctions in the main one.
    """
    pd = job.dir / "phonon"
    pd.mkdir(exist_ok=True)
    ph_scratch = job.scratch / "phonon"
    ph_scratch.mkdir(parents=True, exist_ok=True)
    outdir = str(ph_scratch) + "/"
    atoms = qein.read_structure(text)
    q = [int(x) for x in (opts.get("ph_qgrid") or default_qgrid(atoms))]
    job.log(f"Phonon q-grid: {q[0]}x{q[1]}x{q[2]}")

    def ph_scf():
        t = qein.set_param(text, "control", "calculation", "scf")
        t = qein.set_param(t, "control", "outdir", outdir)
        t = qein.set_param(t, "control", "pseudo_dir", "../pseudos/")
        t = qein.set_param(t, "control", "tprnfor", True)
        cthr = qein.get(qein.parse_input(text), "electrons", "conv_thr", 1e-6)
        try:
            cthr = float(cthr)
        except (TypeError, ValueError):
            cthr = 1e-6
        t = qein.set_param(t, "electrons", "conv_thr", min(cthr, 1e-12 * max(1, len(atoms))))
        (pd / "scf.in").write_text(t)
        rc = mgr.run(job, "pw.x", pd / "scf.in", pd / "scf.out", np_, cwd=pd)
        _check(job, "ph_scf", pd / "scf.out", rc)
    _timed(job, "ph_scf", ph_scf)

    def ph():
        (pd / "ph.in").write_text(phonon.ph_input(prefix, outdir, q))
        rc = mgr.run(job, "ph.x", pd / "ph.in", pd / "ph.out", np_, cwd=pd)
        _check(job, "ph", pd / "ph.out", rc)
    _timed(job, "ph", ph)

    def q2r():
        (pd / "q2r.in").write_text(phonon.q2r_input(prefix))
        rc = mgr.run(job, "q2r.x", pd / "q2r.in", pd / "q2r.out", 1, use_mpi=False, cwd=pd)
        _check(job, "q2r", pd / "q2r.out", rc)
    _timed(job, "q2r", q2r)

    def matdyn():
        kp = qein.make_kpath(atoms, path=(opts.get("kpath") or "").strip() or None, density=25)
        (pd / "qpath.json").write_text(json.dumps(kp, indent=1))
        (pd / "matdyn.in").write_text(phonon.matdyn_path_input(prefix, kp["card"]))
        rc = mgr.run(job, "matdyn.x", pd / "matdyn.in", pd / "matdyn.out", 1, use_mpi=False, cwd=pd)
        _check(job, "matdyn", pd / "matdyn.out", rc)
        grid = qein.kgrid_for_density(atoms, 24 if len(atoms) <= 8 else 16)
        # broadening scaled to the frequency range: a fixed 3 cm^-1 is fine for Si (520 cm^-1) but leaves the
        # DOS of stiff materials (diamond, 1350 cm^-1) spiky on the same mesh
        try:
            wmax = float(np.max(phonon.parse_freq(pd / f"{prefix}.freq")["freq"]))
        except Exception:
            wmax = 0.0
        degauss = round(max(3.0, wmax / 150), 1)
        (pd / "matdyn_dos.in").write_text(phonon.matdyn_dos_input(prefix, grid, degauss_cm=degauss))
        rc = mgr.run(job, "matdyn.x", pd / "matdyn_dos.in", pd / "matdyn_dos.out", 1, use_mpi=False, cwd=pd)
        _check(job, "matdyn", pd / "matdyn_dos.out", rc)
        mesh = qein.kgrid_for_density(atoms, 16 if len(atoms) <= 8 else 8)
        (pd / "matdyn_mesh.in").write_text(phonon.matdyn_mesh_input(prefix, mesh))
        rc = mgr.run(job, "matdyn.x", pd / "matdyn_mesh.in", pd / "matdyn_mesh.out", 1, use_mpi=False, cwd=pd)
        _check(job, "matdyn", pd / "matdyn_mesh.out", rc, required=False)
    _timed(job, "matdyn", matdyn)


def analyse_phonons(d: Path, prefix: str, symbols: list[str]) -> dict | None:
    pd = d / "phonon"
    fq = pd / f"{prefix}.freq"
    if not fq.exists():
        return None
    natoms = len(symbols)
    res: dict = {}
    f = phonon.parse_freq(fq)
    kp = json.loads((pd / "qpath.json").read_text()) if (pd / "qpath.json").exists() else None
    x = qeout.kpath_distance(f["q"], kp["jumps"] if kp else [])
    ticks = []
    if kp:
        for lbl, i in zip(kp["labels"], kp["indices"]):
            if i < len(x):
                name = qein.pretty_label(lbl)
                if ticks and abs(ticks[-1]["x"] - x[i]) < 1e-9:
                    ticks[-1]["label"] += "|" + name
                else:
                    ticks.append({"x": float(x[i]), "label": name})
    res["x"] = x.tolist()
    res["freq"] = f["freq"].T.tolist()            # [branch][q], cm^-1
    res["ticks"] = ticks
    fmin, fmax = float(f["freq"].min()), float(f["freq"].max())
    res["min_cm1"], res["max_cm1"] = fmin, fmax
    # Gamma frequencies after the acoustic sum rule (from the path), and before it (from fildyn1)
    g_idx = next((i for i, q in enumerate(f["q"]) if np.allclose(q, 0, atol=1e-6)), None)
    if g_idx is not None:
        res["gamma_cm1"] = f["freq"][g_idx].tolist()
    dyn1 = pd / f"{prefix}.dyn1"
    if dyn1.exists():
        raw = phonon.parse_dyn_gamma(dyn1)
        res["gamma_raw_cm1"] = raw
        res["asr_violation_cm1"] = float(max(abs(v) for v in sorted(raw, key=abs)[:3])) if len(raw) >= 3 else None
    dos_f = pd / f"{prefix}.phdos"
    mesh_f = pd / f"{prefix}.mesh.freq"
    if mesh_f.exists():
        th = phonon.thermodynamics_modes(phonon.parse_freq(mesh_f)["freq"], natoms)
        res["thermo"] = th
        res["zpe_mev_atom"] = th["zpe_ev_atom"] * 1000
        res["cv_300"] = phonon.value_at(th, "Cv", 300)
        res["s_300"] = phonon.value_at(th, "S", 300)
    if dos_f.exists():
        dd = phonon.parse_phdos(dos_f)
        res["dos"] = {"freq": dd["freq"].tolist(), "total": dd["dos"].tolist(), "elements": {}}
        for i, sym in enumerate(symbols[: len(dd["pdos_atoms"])]):
            cur = res["dos"]["elements"].setdefault(sym, np.zeros_like(dd["dos"]))
            res["dos"]["elements"][sym] = cur + dd["pdos_atoms"][i]
        res["dos"]["elements"] = {k: v.tolist() for k, v in res["dos"]["elements"].items()}
        if "thermo" not in res:  # older runs without the q-mesh: fall back to the DOS
            th = phonon.thermodynamics(dd["freq"], dd["dos"], natoms)
            res["thermo"] = th
            res["zpe_mev_atom"] = th["zpe_ev_atom"] * 1000
            res["cv_300"] = phonon.value_at(th, "Cv", 300)
            res["s_300"] = phonon.value_at(th, "S", 300)
    m = re.search(r"nq1\s*=\s*(\d+),\s*nq2\s*=\s*(\d+),\s*nq3\s*=\s*(\d+)", (pd / "ph.in").read_text()) if (pd / "ph.in").exists() else None
    res["qgrid"] = [int(m.group(i)) for i in (1, 2, 3)] if m else None
    res["stable"] = fmin > phonon.IMAG_THRESHOLD_CM1
    return res


# --------------------------------------------------------------------------- IR / Raman
def run_spectra(mgr: Manager, job: Job, text: str, prefix: str, scf_res: dict, np_: int):
    """Gamma-point ph.x with epsil (and lraman when allowed), then dynmat.x for TO and LO modes."""
    sd = job.dir / "spectra"
    sd.mkdir(exist_ok=True)
    parsed = qein.parse_input(text)
    occ = str(qein.get(parsed, "system", "occupations", "fixed")).lower()
    if occ != "fixed":
        msg = "IR/Raman need an insulator with occupations = 'fixed' (the dielectric response of a metal is not defined)."
        for s_ in ("sp_scf", "sp_ph", "dynmat"):
            job.step(s_, status="skipped", error=msg)
        (sd / "skipped.txt").write_text(msg)
        return
    types = sorted({cif.pseudo_kind(p) for p in (job.dir / "pseudos").iterdir() if p.is_file()})
    raman, why = spectra.raman_allowed(scf_res.get("xc"), types, occ)
    (sd / "raman_allowed.json").write_text(json.dumps({"raman": raman, "reason": why}))
    if not raman:
        job.log(f"Raman skipped: {why} IR is still calculated.")
    scratch = job.scratch / "spectra"
    scratch.mkdir(parents=True, exist_ok=True)
    outdir = str(scratch) + "/"

    def sp_scf():
        t = qein.set_param(text, "control", "calculation", "scf")
        t = qein.set_param(t, "control", "outdir", outdir)
        t = qein.set_param(t, "control", "pseudo_dir", "../pseudos/")
        t = qein.set_param(t, "electrons", "conv_thr", 1e-12)
        (sd / "scf.in").write_text(t)
        rc = mgr.run(job, "pw.x", sd / "scf.in", sd / "scf.out", np_, cwd=sd)
        _check(job, "sp_scf", sd / "scf.out", rc)
    _timed(job, "sp_scf", sp_scf)

    def sp_ph():
        (sd / "ph.in").write_text(spectra.ph_gamma_input(prefix, outdir, raman))
        rc = mgr.run(job, "ph.x", sd / "ph.in", sd / "ph.out", np_, cwd=sd)
        _check(job, "sp_ph", sd / "ph.out", rc)
    _timed(job, "sp_ph", sp_ph)

    def dyn():
        (sd / "dynmat.in").write_text(spectra.dynmat_input(prefix))
        rc = mgr.run(job, "dynmat.x", sd / "dynmat.in", sd / "dynmat.out", 1, use_mpi=False, cwd=sd)
        _check(job, "dynmat", sd / "dynmat.out", rc)
        (sd / "dynmat_lo.in").write_text(spectra.dynmat_input(prefix, [1.0, 0.0, 0.0]))
        rc = mgr.run(job, "dynmat.x", sd / "dynmat_lo.in", sd / "dynmat_lo.out", 1, use_mpi=False, cwd=sd)
        _check(job, "dynmat", sd / "dynmat_lo.out", rc, required=False)
    _timed(job, "dynmat", dyn)


# --------------------------------------------------------------------------- UV-Vis
def run_uvvis(mgr: Manager, job: Job, text: str, prefix: str, opts: dict, np_: int):
    """SCF + NSCF with many empty bands (symmetry-reduced k-grid), then the dielectric function from the
    HDF5 wavefunctions (qestudio/optics.py)."""
    ud = job.dir / "uvvis"
    ud.mkdir(exist_ok=True)
    types = sorted({cif.pseudo_kind(p) for p in (job.dir / "pseudos").iterdir() if p.is_file()})
    if any(t not in ("NC", "SL") for t in types):
        msg = ("UV-Vis needs norm-conserving pseudopotentials (ultrasoft/PAW need augmentation terms in the optical "
               "matrix elements). Convert the CIF with 'LDA norm-conserving' or use NC pseudopotentials.")
        for s_ in ("uv_nscf", "optics"):
            job.step(s_, status="skipped", error=msg)
        (ud / "skipped.txt").write_text(msg)
        return
    scratch = job.scratch / "uvvis"
    scratch.mkdir(parents=True, exist_ok=True)
    outdir = str(scratch) + "/"
    atoms = qein.read_structure(text)
    parsed = qein.parse_input(text)
    noncolin = bool(qein.get(parsed, "system", "noncolin", False))
    scf_out = qeout.parse_pw((job.dir / "scf.out").read_text(errors="replace"))
    nelec = scf_out.get("nelec") or 0
    nocc = int(math.ceil(nelec if noncolin else nelec / 2))
    nbnd = int(opts.get("uv_nbnd") or (2 * nocc + 8))
    grid = [int(x) for x in (opts.get("uv_grid") or qein.kgrid_for_density(atoms, 16 if len(atoms) <= 4 else 10))]

    def uv_nscf():
        t = qein.set_param(text, "control", "calculation", "scf")
        t = qein.set_param(t, "control", "outdir", outdir)
        t = qein.set_param(t, "control", "pseudo_dir", "../pseudos/")
        (ud / "scf.in").write_text(t)
        rc = mgr.run(job, "pw.x", ud / "scf.in", ud / "scf.out", np_, cwd=ud)
        _check(job, "uv_nscf", ud / "scf.out", rc)
        t = qein.set_param(t, "control", "calculation", "nscf")
        t = qein.set_param(t, "system", "nbnd", nbnd)
        t = qein.replace_card(t, "K_POINTS", "K_POINTS automatic\n  {} {} {} 0 0 0\n".format(*grid))
        (ud / "nscf.in").write_text(t)
        rc = mgr.run(job, "pw.x", ud / "nscf.in", ud / "nscf.out", np_, cwd=ud)
        _check(job, "uv_nscf", ud / "nscf.out", rc)
        job.log(f"UV-Vis NSCF: {grid[0]}x{grid[1]}x{grid[2]} k-grid, {nbnd} bands")
    _timed(job, "uv_nscf", uv_nscf)

    def opt():
        save = scratch / f"{prefix}.save"
        if not list(save.glob("wfc*.hdf5")):
            raise RuntimeError("No HDF5 wavefunctions found (QE must be built with HDF5, as the conda-forge build is).")
        meta = optics.compute(save, ud, emax_ev=float(opts.get("uv_emax") or 20), sigma_ev=float(opts.get("uv_sigma") or 0.1),
                              shift_ev=float(opts.get("uv_scissor") or 0))
        meta.update(grid=grid, pseudo_types=types)
        (ud / "optics.json").write_text(json.dumps(meta, indent=1))
        if meta.get("metal"):
            job.log(f"UV-Vis (metal): Drude plasma energy {meta['drude_plasma_ev']:.2f} eV, f-sum ratio {meta['fsum_ratio']:.2f}")
        else:
            job.log(f"UV-Vis: eps1(0) = {meta['eps_static']:.2f}, f-sum ratio {meta['fsum_ratio']:.2f}")
    _timed(job, "optics", opt)


def analyse_uvvis(d: Path) -> dict | None:
    ud = d / "uvvis"
    if (ud / "skipped.txt").exists():
        return {"skipped": (ud / "skipped.txt").read_text()}
    if not (ud / "optics.dat").exists():
        return None
    data = np.loadtxt(ud / "optics.dat")
    meta = json.loads((ud / "optics.json").read_text())
    E, e1, e2, n, k, alpha, R = data.T
    step = max(1, len(E) // 1500)                          # keep the payload small
    res = {"meta": meta, "E": E[::step].tolist(), "eps1": e1[::step].tolist(), "eps2": e2[::step].tolist(),
           "n": n[::step].tolist(), "k": k[::step].tolist(), "alpha": alpha[::step].tolist(), "R": R[::step].tolist()}
    i589 = int(np.argmin(abs(E - optics.HC_EV_NM / 589.0)))
    res["n_589"] = float(n[i589])
    res["k_589"] = float(k[i589])
    res["R_589"] = float(R[i589])
    if meta.get("metal"):
        # a metal absorbs at every energy: no onset; the screened plasma edge is where eps1 crosses zero upwards
        up = np.nonzero((e1[:-1] < 0) & (e1[1:] >= 0) & (E[:-1] > 1.0))[0]
        if len(up):
            res["plasma_edge_ev"] = float(E[up[-1] + 1])
        inter = E > 0.3          # the interband peak, not the Drude divergence at w -> 0
        res["peak_eps2_ev"] = float(E[inter][int(np.argmax(e2[inter] - optics.drude(E[inter], meta["drude_plasma_ev"],
                                                                                         meta["drude_gamma_ev"])[1]))])
        return res
    on = np.nonzero(alpha >= 1e4)[0]
    if len(on):
        res["onset_ev"] = float(E[on[0]])
        res["onset_nm"] = float(optics.HC_EV_NM / E[on[0]])
    res["peak_eps2_ev"] = float(E[int(np.argmax(e2))])
    return res


# --------------------------------------------------------------------------- NMR
def build_gipaw(job: Job, out: Path):
    """Run build_gipaw.sh (compilers, QE source, pw.x + gipaw.x into ~/gipaw/bin). Cancel stops it."""
    script = Path(__file__).resolve().parent.parent / "build_gipaw.sh"
    job.log(f"GIPAW is not installed yet: building it once ($ bash {script} > {out.name})")
    with open(out, "w") as fo:
        code = script.read_text().replace("\r\n", "\n")       # a Windows checkout may have CRLF endings
        job.proc = subprocess.Popen(["bash", "-c", code, str(script)], stdout=fo, stderr=subprocess.STDOUT,
                                    stdin=subprocess.DEVNULL, start_new_session=True)
        rc = job.proc.wait()
    if job.cancelled:
        raise RuntimeError("Cancelled by user")
    if rc != 0 or not nmr.available():
        tail = out.read_text(errors="replace").strip().splitlines()[-3:]
        job.step("gipaw_build", status="failed", error="GIPAW build failed: " + " | ".join(tail))
        job.log("GIPAW build failed: " + " | ".join(tail))


def run_nmr(mgr: Manager, job: Job, text: str, prefix: str, np_: int):
    """SCF and gipaw.x (job = 'nmr') with the self-built QE-GIPAW, in runs/<id>/nmr/."""
    nd = job.dir / "nmr"
    nd.mkdir(exist_ok=True)

    has_build_step = any(s.get("name") == "gipaw_build" for s in job.load().get("steps", []))

    def skip(msg):
        done = {s["name"] for s in job.load().get("steps", []) if s.get("status") in ("done", "failed")}
        for s_ in ("gipaw_build", "nmr_scf", "nmr"):
            if s_ not in done:
                job.step(s_, status="skipped", error=msg)
        (nd / "skipped.txt").write_text(msg)

    # cheap checks first, so a 30-60 minute build is never started for a run that cannot use it
    occ = str(qein.get(qein.parse_input(text), "system", "occupations", "fixed")).lower()
    if occ != "fixed":
        return skip("NMR shielding (GIPAW job 'nmr') needs an insulator with occupations = 'fixed'.")
    bad = [f"{p.name} ({why})" for p in sorted((job.dir / "pseudos").iterdir())
           if p.is_file() for ok, why in [nmr.gipaw_ready(p)] if not ok]
    if bad:
        return skip("These pseudopotentials have no usable GIPAW data: " + "; ".join(bad) +
                    ". Use PSlibrary PAW files (CIF converter: 'PBE PAW, for NMR').")
    if not nmr.available():
        if config.PUBLIC or not has_build_step:
            return skip("GIPAW is not installed on this computer. Run build_gipaw.sh once in Ubuntu "
                        "(README, installation step 5; about 30-60 minutes, no admin rights needed), then rerun.")
        _timed(job, "gipaw_build", lambda: build_gipaw(job, nd / "gipaw_build.out"))
        if not nmr.available():
            return skip("The automatic GIPAW build did not finish; see nmr/gipaw_build.out and ~/gipaw/*.log.")
    elif has_build_step:
        job.step("gipaw_build", status="done", seconds=0)   # built meanwhile (e.g. by an earlier job)
    scratch = job.scratch / "nmr"
    scratch.mkdir(parents=True, exist_ok=True)
    outdir = str(scratch) + "/"

    def nmr_scf():
        t = qein.set_param(text, "control", "calculation", "scf")
        t = qein.set_param(t, "control", "outdir", outdir)
        t = qein.set_param(t, "control", "pseudo_dir", "../pseudos/")
        t = qein.set_param(t, "electrons", "conv_thr", 1e-10)
        (nd / "scf.in").write_text(t)
        rc = mgr.run(job, "pw.x", nd / "scf.in", nd / "scf.out", np_, cwd=nd, gipaw=True)
        _check(job, "nmr_scf", nd / "scf.out", rc)
    _timed(job, "nmr_scf", nmr_scf)

    def gip():
        (nd / "gipaw.in").write_text(nmr.gipaw_input(prefix, outdir))
        rc = mgr.run(job, "gipaw.x", nd / "gipaw.in", nd / "gipaw.out", np_, cwd=nd, gipaw=True)
        out = (nd / "gipaw.out").read_text(errors="replace") if (nd / "gipaw.out").exists() else ""
        if "JOB DONE" not in out and "S matrix not positive definite" in out:
            # Davidson's subspace rotation at k+q breaks down for some PAW sets with semicore states
            # (NaF with Na.pbe-spn-kjpaw); conjugate gradients works, several times slower
            job.log("gipaw.x: 'S matrix not positive definite' with Davidson; retrying with diagonalization = 'cg'")
            (nd / "gipaw_david.out").write_text(out)
            (nd / "gipaw.in").write_text(nmr.gipaw_input(prefix, outdir, diagonalization="cg"))
            rc = mgr.run(job, "gipaw.x", nd / "gipaw.in", nd / "gipaw.out", np_, cwd=nd, gipaw=True)
        _check(job, "nmr", nd / "gipaw.out", rc)
    _timed(job, "nmr", gip)


# --------------------------------------------------------------------------- convergence
def convergence(mgr: Manager, job: Job, text: str, kind: str, opts: dict, np_: int):
    parsed = qein.parse_input(text)
    ecut0 = float(qein.get(parsed, "system", "ecutwfc"))
    rho0 = qein.get(parsed, "system", "ecutrho")
    ratio = float(rho0) / ecut0 if rho0 else 4.0
    kinfo = qein.kpoints_info(parsed)
    atoms = qein.read_structure(text)
    grids: dict = {}
    if kind == "ecut":
        vals = [float(v) for v in re.split(r"[,\s]+", str(opts.get("ecut_list") or "")) if v] or \
            [round(ecut0 * f) for f in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)]
    else:
        if not kinfo.get("grid"):
            job.step("conv_k", status="skipped", error="K_POINTS is not 'automatic'")
            return
        user = [int(v) for v in re.split(r"[,\s]+", str(opts.get("k_list") or "")) if v]
        if user:      # "k (largest)" values: grids of that density along the longest reciprocal vector
            for v in user:
                grids[v] = qein.kgrid_for_density(atoms, v)
        else:         # scale the input's own grid, so the 1x point is exactly the grid being used
            g0 = [int(x) for x in kinfo["grid"]]
            for f in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5):
                g = [1 if x == 1 else max(1, int(round(x * f))) for x in g0]   # 1 stays 1 (slab / vacuum direction)
                if max(g) not in grids:               # skip duplicates from rounding small grids
                    grids[max(g)] = g
        vals = sorted(grids)
    root = job.dir / f"conv_{kind}"
    root.mkdir(exist_ok=True)
    rows = []
    for v in vals:
        if job.cancelled:
            raise RuntimeError("Cancelled by user")
        sub = root / (f"ecut_{v:g}" if kind == "ecut" else f"k_{v}")
        sub.mkdir(exist_ok=True)
        t = qein.set_param(text, "control", "calculation", "scf")
        sub_scratch = job.scratch / f"conv_{kind}" / sub.name
        sub_scratch.mkdir(parents=True, exist_ok=True)  # pw.x only creates the last directory level
        t = qein.set_param(t, "control", "outdir", str(sub_scratch) + "/")
        t = qein.set_param(t, "control", "pseudo_dir", "../../pseudos/")
        t = qein.set_param(t, "control", "tstress", True)
        t = qein.set_param(t, "control", "tprnfor", True)
        t = qein.set_param(t, "control", "disk_io", "nowf")
        if kind == "ecut":
            t = qein.set_param(t, "system", "ecutwfc", float(v))
            t = qein.set_param(t, "system", "ecutrho", float(v) * ratio)
            label = f"ecutwfc = {v:g} Ry"
        else:
            grid = grids[v]
            sh = kinfo.get("shift", [0, 0, 0])
            t = qein.replace_card(t, "K_POINTS", "K_POINTS automatic\n  {} {} {} {} {} {}\n".format(*grid, *sh))
            label = "k = {}x{}x{}".format(*grid)
        (sub / "scf.in").write_text(t)
        job.log(f"convergence {kind}: {label}")
        rc = mgr.run(job, "pw.x", sub / "scf.in", sub / "scf.out", np_, cwd=sub)
        out = (sub / "scf.out").read_text(errors="replace")
        r = qeout.parse_pw(out)
        shutil.rmtree(sub_scratch, ignore_errors=True)
        rows.append({"value": v, "label": label, "energy_ry": r.get("energy_ry"), "pressure_kbar": r.get("pressure_kbar"),
                     "total_force": r.get("total_force"), "wall_s": r.get("wall_s"), "nks": r.get("nks"),
                     "ok": r.get("job_done") and not r.get("not_converged"),
                     "grid": grid if kind == "k" else None})
    (root / "results.json").write_text(json.dumps(rows, indent=1))
    bad = [r["label"] for r in rows if not r["ok"]]
    if bad:
        first = next((s for s in sorted(root.glob("*/scf.out")) if "JOB DONE" not in s.read_text(errors="replace")), None)
        err = "; ".join(qeout.errors_of(first.read_text(errors="replace"))) if first else ""
        job.step(f"conv_{kind}", status="failed",
                 error=f"{len(bad)}/{len(rows)} runs failed ({', '.join(bad)}). {err}".strip())


def convergence_summary(rows: list[dict], nat: int, thr_mev: float = 1.0) -> dict:
    good = [r for r in rows if r.get("energy_ry") is not None]
    if len(good) < 2:
        return {"rows": rows}
    ref = good[-1]["energy_ry"]
    for r in good:
        r["delta_mev_atom"] = (r["energy_ry"] - ref) * qeout.RY_EV * 1000 / max(nat, 1)
        if r.get("pressure_kbar") is not None and good[-1].get("pressure_kbar") is not None:
            r["delta_p_kbar"] = r["pressure_kbar"] - good[-1]["pressure_kbar"]
    rec = None
    for i, r in enumerate(good[:-1]):
        if all(abs(x["delta_mev_atom"]) < thr_mev for x in good[i:]):
            rec = r
            break
    return {"rows": rows, "threshold_mev_atom": thr_mev, "recommended": rec["value"] if rec else None,
            "recommended_label": rec["label"] if rec else None}


# --------------------------------------------------------------------------- analysis
def analyse_job(job: Job) -> dict:
    d = job.dir
    res: dict = {"id": job.id}
    st = job.load()
    scf_text = (d / "scf.in").read_text()
    parsed = qein.parse_input(scf_text)
    occ = str(qein.get(parsed, "system", "occupations", "fixed")).lower()
    noncolin = bool(qein.get(parsed, "system", "noncolin", False))
    try:
        atoms = qein.read_structure(scf_text)
        res["structure"] = {"formula": atoms.get_chemical_formula(mode="hill"), "natoms": len(atoms),
                            "symbols": atoms.get_chemical_symbols(), "positions": atoms.get_positions().tolist(),
                            "cell": atoms.cell.array.tolist(), "abc": list(atoms.cell.cellpar()),
                            "volume_A3": atoms.get_volume()}
    except Exception as e:
        res["structure_error"] = str(e)
    scf = qeout.parse_pw((d / "scf.out").read_text(errors="replace")) if (d / "scf.out").exists() else {}
    res["scf"] = _slim(scf)
    nat = scf.get("nat") or res.get("structure", {}).get("natoms") or 1
    fixed = occ != "smearing"
    summary = {
        "qe_version": scf.get("version"), "calculation": scf.get("calculation"),
        "energy_ry": scf.get("energy_ry"),
        "energy_ev": scf["energy_ry"] * qeout.RY_EV if scf.get("energy_ry") is not None else None,
        "energy_ev_atom": scf["energy_ry"] * qeout.RY_EV / nat if scf.get("energy_ry") is not None else None,
        "fermi_ev": scf.get("fermi_ev"), "homo_ev": scf.get("homo_ev"), "lumo_ev": scf.get("lumo_ev"),
        "scf_iterations": scf.get("scf_iterations"), "wall_s": scf.get("wall_s"),
        "total_force": scf.get("total_force"), "pressure_kbar": scf.get("pressure_kbar"),
        "total_mag": scf.get("total_mag"), "abs_mag": scf.get("abs_mag"), "nelec": scf.get("nelec"),
        "nks": scf.get("nks"), "ecutwfc": scf.get("ecutwfc"), "ecutrho": scf.get("ecutrho"),
        "xc": scf.get("xc"), "mpi_procs": scf.get("mpi_procs"), "occupations": occ,
    }
    # gap from the SCF grid
    sb = qeout.bands_from_pw(scf)
    ref_scf = scf.get("fermi_ev")
    if sb:
        summary["gap_scf_grid"] = _jsonable(qeout.band_edges(sb["eig"], sb["kpts"], scf.get("nelec"), ref_scf, fixed, noncolin))

    # band structure
    if (d / "bands.out").exists():
        bo = qeout.parse_pw((d / "bands.out").read_text(errors="replace"))
        bd = qeout.bands_from_pw(bo)
        if bd:
            kp = json.loads((d / "kpath.json").read_text()) if (d / "kpath.json").exists() else None
            jumps = kp["jumps"] if kp else []
            x = qeout.kpath_distance(bd["kpts"], jumps)
            edges = qeout.band_edges(bd["eig"], bd["kpts"], scf.get("nelec"), ref_scf, fixed, noncolin)
            ticks = []
            if kp:
                for lbl, i in zip(kp["labels"], kp["indices"]):
                    if i < len(x):
                        name = qein.pretty_label(lbl)
                        if ticks and abs(ticks[-1]["x"] - x[i]) < 1e-9:
                            ticks[-1]["label"] += "|" + name
                        else:
                            ticks.append({"x": float(x[i]), "label": name})
            else:
                ticks = [{"x": float(x[i]), "label": ""} for i in qeout.detect_kinks(bd["kpts"])]
            # describe where VBM/CBM are
            if edges.get("metal") is False:
                for which in ("vbm", "cbm"):
                    i = edges[f"{which}_index"]
                    edges[f"{which}_x"] = float(x[i])
                    edges[f"{which}_where"] = _where(x[i], ticks)
            res["bands"] = {"x": x.tolist(), "kpts": bd["kpts"].tolist(),
                            "eig": [e.T.tolist() for e in bd["eig"]],  # per spin: [band][k]
                            "ticks": ticks, "edges": _jsonable(edges), "nbnd": bd["eig"][0].shape[1]}
            summary["gap_bands"] = _jsonable(edges)
    # nscf (dense grid) gap
    if (d / "nscf.out").exists():
        no = qeout.parse_pw((d / "nscf.out").read_text(errors="replace"))
        nb = qeout.bands_from_pw(no)
        if nb:
            summary["gap_nscf"] = _jsonable(qeout.band_edges(nb["eig"], nb["kpts"], scf.get("nelec"),
                                                             no.get("fermi_ev", ref_scf), fixed, noncolin))
            summary["fermi_nscf_ev"] = no.get("fermi_ev")
            summary["homo_nscf_ev"] = no.get("homo_ev")
    # best overall gap: smallest over all k-sets
    # only sources that could decide metal/insulator (e.g. the SCF grid has no conduction band when nbnd = nocc)
    gaps = [summary.get(k) for k in ("gap_bands", "gap_nscf", "gap_scf_grid")
            if summary.get(k) and summary[k].get("metal") is not None]
    if gaps:
        if any(g.get("metal") for g in gaps):
            summary["is_metal"] = True
            summary["gap_ev"] = 0.0
        else:
            vbm = max(g["vbm_ev"] for g in gaps if g.get("vbm_ev") is not None)
            cbm = min(g["cbm_ev"] for g in gaps if g.get("cbm_ev") is not None)
            summary["is_metal"] = False
            summary["gap_ev"] = cbm - vbm
            summary["vbm_ev"], summary["cbm_ev"] = vbm, cbm
            if summary.get("gap_bands"):
                summary["gap_type"] = summary["gap_bands"].get("gap_type")
                summary["direct_gap_ev"] = summary["gap_bands"].get("direct_gap_ev")
    # energy reference for plots: VBM for insulators, E_F for metals
    if summary.get("is_metal"):
        ref = summary.get("fermi_nscf_ev") or ref_scf or (summary.get("gap_bands") or {}).get("reference_ev")
    else:
        ref = summary.get("vbm_ev") if summary.get("vbm_ev") is not None else (scf.get("homo_ev") or ref_scf)
    summary["energy_reference_ev"] = ref
    summary["energy_reference_label"] = "E_F" if summary.get("is_metal") else "VBM"

    # DOS
    prefix = st.get("prefix", "pwscf")
    if (d / f"{prefix}.dos").exists():
        try:
            res["dos"] = qeout.parse_dos_file(d / f"{prefix}.dos")
        except Exception as e:
            res["dos_error"] = str(e)
    if list(d.glob(f"{prefix}.pdos_atm#*")):
        try:
            res["pdos"] = qeout.parse_pdos(d, prefix)
        except Exception as e:
            res["pdos_error"] = str(e)
    if (d / "projwfc.out").exists():
        res["projwfc"] = qeout.parse_projwfc_out((d / "projwfc.out").read_text(errors="replace"))
    # bands.x gnu file
    gnu = d / f"{prefix}_bands.dat.gnu"
    if gnu.exists():
        try:
            data = np.loadtxt(gnu)
            res["bandsx_gnu"] = {"file": gnu.name, "rows": int(data.shape[0])}
        except Exception:
            pass
    # convergence
    for kind in ("ecut", "k"):
        f = d / f"conv_{kind}" / "results.json"
        if f.exists():
            res[f"conv_{kind}"] = convergence_summary(json.loads(f.read_text()), nat,
                                                      float(st["options"].get("conv_thr_mev") or 1.0))
    # phonons
    try:
        ph = analyse_phonons(d, prefix, res.get("structure", {}).get("symbols", []))
        if ph:
            res["phonons"] = ph
    except Exception as e:
        res["phonons_error"] = str(e)
    # NMR
    nd = d / "nmr"
    if (nd / "skipped.txt").exists():
        res["nmr"] = {"skipped": (nd / "skipped.txt").read_text()}
    elif (nd / "gipaw.out").exists():
        try:
            res["nmr"] = nmr.parse_gipaw((nd / "gipaw.out").read_text(errors="replace"))
        except Exception as e:
            res["nmr_error"] = str(e)
    # UV-Vis
    try:
        uv = analyse_uvvis(d)
        if uv is not None:
            res["uvvis"] = uv
    except Exception as e:
        res["uvvis_error"] = str(e)
    # IR / Raman
    try:
        sp = spectra.analyse(d, prefix)
        if sp is not None:
            ra = d / "spectra" / "raman_allowed.json"
            if ra.exists():
                sp.update(json.loads(ra.read_text()))
            res["spectra"] = sp
        elif (d / "spectra" / "skipped.txt").exists():
            res["spectra"] = {"skipped": (d / "spectra" / "skipped.txt").read_text()}
    except Exception as e:
        res["spectra_error"] = str(e)
    # checks
    checks = []
    for out in (sorted(d.glob("*.out")) + sorted(d.glob("conv_*/*/scf.out")) + sorted(d.glob("phonon/*.out"))
                + sorted(d.glob("spectra/*.out")) + sorted(d.glob("uvvis/*.out")) + sorted(d.glob("nmr/*.out"))):
        txt = out.read_text(errors="replace")
        rel = str(out.relative_to(d))
        if "JOB DONE" not in txt:
            checks.append({"level": "error", "msg": f"{rel} does not end with JOB DONE."})
        if "convergence NOT achieved" in txt:
            checks.append({"level": "error", "msg": f"{rel}: convergence NOT achieved."})
    if summary.get("total_force") and summary["total_force"] > 1e-3 and scf.get("calculation") == "scf":
        checks.append({"level": "warn", "msg": f"Total force {summary['total_force']:.4f} Ry/bohr: the structure is not relaxed."})
    if summary.get("pressure_kbar") is not None and abs(summary["pressure_kbar"]) > 10:
        checks.append({"level": "info", "msg": f"Pressure {summary['pressure_kbar']:.1f} kbar: the lattice constant is not at equilibrium for this functional/cutoff."})
    if (d / "relax.json").exists():
        rx = json.loads((d / "relax.json").read_text())
        res["relax"] = rx
        ab, aa = rx.get("abc_before") or [], rx.get("abc_after") or []
        dv = (rx["volume_after_A3"] / rx["volume_before_A3"] - 1) * 100 if rx.get("volume_before_A3") else 0
        lens = ", ".join(f"{b:.4f} → {a:.4f}" for b, a in zip(ab, aa)) if rx["mode"] == "vc-relax" else ""
        fp = (f"force {rx['force_first']:.4f} → {rx['force']:.5f} Ry/bohr" if rx.get("force") is not None
              and rx.get("force_first") is not None else "")
        pp = (f", pressure {rx['pressure_kbar_first']:.1f} → {rx['pressure_kbar']:.1f} kbar"
              if rx["mode"] == "vc-relax" and rx.get("pressure_kbar") is not None else "")
        msg = (f"Structure relaxed first ({rx['mode']}, {rx['bfgs_steps']} SCF cycles): "
               + (f"a, b, c (Å) {lens}; volume {dv:+.1f}%; " if lens else "") + fp + pp
               + ". All results below are for the relaxed structure (scf.in).")
        checks.append({"level": "info" if rx.get("converged") else "warn",
                       "msg": msg if rx.get("converged") else msg + " The relaxation did NOT converge: the last structure was used."})
    if res.get("projwfc", {}).get("spilling") and res["projwfc"]["spilling"] > 0.05:
        checks.append({"level": "info", "msg": f"PDOS spilling parameter {res['projwfc']['spilling']:.3f}: projections are approximate."})
    # metals: a Fermi level that moves on the denser NSCF grid means the SCF k-grid is not converged
    ef_scf, ef_nscf = summary.get("fermi_ev"), summary.get("fermi_nscf_ev")
    if summary.get("is_metal") and ef_scf is not None and ef_nscf is not None and abs(ef_nscf - ef_scf) > 0.05:
        checks.append({"level": "warn", "msg": (
            f"Fermi energy changes from {ef_scf:.3f} eV (SCF grid) to {ef_nscf:.3f} eV on the denser NSCF grid "
            f"({(ef_nscf - ef_scf) * 1000:+.0f} meV): the SCF k-point grid is not converged for this metal. "
            "Plots use the NSCF value. Use a denser k-grid and/or a larger degauss, and run the k-point convergence test.")})
    for kind, name in (("ecut", "cutoff"), ("k", "k-point")):
        conv = res.get(f"conv_{kind}")
        if conv and len([r for r in conv.get("rows", []) if r.get("energy_ry") is not None]) >= 2 and conv.get("recommended") is None:
            checks.append({"level": "warn", "msg": (
                f"The {name} convergence test did not converge to {conv.get('threshold_mev_atom', 1)} meV/atom within the "
                "tested values. Test larger values before trusting total energies.")})
    ph = res.get("phonons")
    if ph:
        if not ph["stable"]:
            checks.append({"level": "warn", "msg": (
                f"Imaginary phonon frequencies (down to {ph['min_cm1']:.0f} cm⁻¹): the structure is dynamically unstable, "
                "not relaxed, or the q-grid / cutoff / conv_thr are not converged. Thermodynamics ignore these modes.")})
        if ph.get("asr_violation_cm1") and ph["asr_violation_cm1"] > 30:
            checks.append({"level": "info", "msg": (
                f"Acoustic modes at Γ were {ph['asr_violation_cm1']:.0f} cm⁻¹ before the acoustic sum rule was applied: "
                "consider a higher ecutwfc or a denser k-grid for more accurate phonons.")})
        if summary.get("total_force") and summary["total_force"] > 1e-3:
            checks.append({"level": "warn", "msg": "Forces are not zero: phonons should be calculated for a relaxed structure."})
    sp = res.get("spectra") or {}
    if sp.get("born_sum") is not None and abs(sp["born_sum"]) > 0.1:
        checks.append({"level": "info", "msg": (
            f"Born effective charges sum to {sp['born_sum']:.3f} instead of 0 (acoustic sum rule); the IR intensities apply "
            "the sum rule, but a higher ecutwfc or denser k-grid would make the raw charges more accurate.")})
    if sp.get("reason"):
        checks.append({"level": "info", "msg": f"Raman not calculated: {sp['reason']} IR was calculated."})
    if sp.get("skipped"):
        checks.append({"level": "info", "msg": sp["skipped"]})
    nm = res.get("nmr") or {}
    if nm.get("skipped"):
        checks.append({"level": "info", "msg": nm["skipped"]})
    uv = res.get("uvvis") or {}
    if uv.get("skipped"):
        checks.append({"level": "info", "msg": uv["skipped"]})
    elif uv.get("meta") and uv["meta"].get("fsum_ratio") is not None and uv["meta"]["fsum_ratio"] < 0.8:
        checks.append({"level": "warn", "msg": (
            f"UV-Vis f-sum rule only {uv['meta']['fsum_ratio']*100:.0f}% fulfilled: too few empty bands for energies up to "
            f"{uv['meta']['emax_ev']:.0f} eV. The low-energy part is still fine.")})
    res["checks"] = checks
    res["summary"] = _jsonable(summary)
    tmp = d / "results.json.tmp"
    tmp.write_text(json.dumps(_jsonable(res)))
    tmp.replace(d / "results.json")  # atomic: the UI may be polling
    _write_summary_csv(d, summary)
    return res


def _where(x: float, ticks: list[dict]) -> str:
    for a, b in zip(ticks, ticks[1:]):
        if a["x"] - 1e-9 <= x <= b["x"] + 1e-9:
            if abs(x - a["x"]) < 1e-6:
                return f"at {a['label']}"
            if abs(x - b["x"]) < 1e-6:
                return f"at {b['label']}"
            frac = (x - a["x"]) / (b["x"] - a["x"]) if b["x"] > a["x"] else 0
            return f"{frac*100:.0f}% of the way {a['label']}→{b['label']}"
    return ""


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


def _write_summary_csv(d: Path, summary: dict):
    import csv
    flat = {}
    for k, v in summary.items():
        if isinstance(v, dict):
            for kk, vv in v.items():
                if not isinstance(vv, (list, dict)):
                    flat[f"{k}.{kk}"] = vv
        else:
            flat[k] = v
    with open(d / "summary.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "value"])
        for k, v in flat.items():
            w.writerow([k, v])


# --------------------------------------------------------------------------- listing
def list_jobs() -> list[dict]:
    out = []
    if not RUNS.exists():
        return out
    for d in sorted(RUNS.iterdir(), reverse=True):
        st = Job(d.name).load()
        if not st:
            continue
        item = {k: st.get(k) for k in ("id", "prefix", "input_name", "status", "created", "finished", "current_step")}
        r = d / "results.json"
        if r.exists():
            try:
                full = json.loads(r.read_text())
                s = full.get("summary", {})
                item.update(energy_ry=s.get("energy_ry"), gap_ev=s.get("gap_ev"), is_metal=s.get("is_metal"))
                if (full.get("nmr") or {}).get("mean_sigma"):
                    item["nmr_sigma"] = full["nmr"]["mean_sigma"]
            except json.JSONDecodeError:
                pass
        out.append(item)
    return out
