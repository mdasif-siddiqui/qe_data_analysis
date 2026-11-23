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

from . import config, pseudos, qein, qeout

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
            cwd: Path | None = None) -> int:
        """Run a QE program, input via -in (pw.x) or stdin (post-processing)."""
        cwd = cwd or job.dir
        env = {**os.environ, "OMP_NUM_THREADS": "1",
               "PATH": f"{QE_BIN}:{os.environ.get('PATH', '')}"}
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
                launcher = [exe("mpirun"), "-np", str(np_)]
                if os.geteuid() == 0:
                    launcher.insert(1, "--allow-run-as-root")
            if mode == "mpi_cma":
                env_m["OMPI_MCA_smsc"] = "^cma"
                env_m["OMPI_MCA_btl_vader_single_copy_mechanism"] = "none"
            cmd = launcher + [exe(program)]
            stdin = subprocess.DEVNULL
            if program == "pw.x":
                cmd += ["-in", inp.name]
            else:
                stdin = open(inp)
            job.log(f"$ {'OMP_NUM_THREADS=1 ' + ' '.join(cmd)}{'' if program == 'pw.x' else ' < ' + inp.name} > {out.name}")
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
    steps = [{"name": "prepare", "label": "Prepare inputs & pseudopotentials"},
             {"name": "scf", "label": "SCF (pw.x)"}]
    if options.get("bands", True):
        steps += [{"name": "bands", "label": "Band structure (pw.x bands)"},
                  {"name": "bandsx", "label": "bands.x post-processing"}]
    if options.get("dos", True):
        steps += [{"name": "nscf", "label": "NSCF on dense k-grid (pw.x)"},
                  {"name": "dos", "label": "Density of states (dos.x)"},
                  {"name": "pdos", "label": "Projected DOS (projwfc.x)"}]
    if options.get("conv"):
        steps += [{"name": "conv_ecut", "label": "Cutoff convergence test"},
                  {"name": "conv_k", "label": "k-point convergence test"}]
    steps.append({"name": "analysis", "label": "Analysis"})
    for s in steps:
        s["status"] = "pending"
    job.save({"id": jid, "prefix": prefix, "input_name": scf_name, "status": "queued",
              "created": time.time(), "options": options, "steps": steps, "notes": notes})
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
        ("S matrix not positive definite", "atoms may be too close: check positions / cell."),
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
        if calc not in ("scf", "relax", "vc-relax"):
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


# --------------------------------------------------------------------------- convergence
def convergence(mgr: Manager, job: Job, text: str, kind: str, opts: dict, np_: int):
    parsed = qein.parse_input(text)
    ecut0 = float(qein.get(parsed, "system", "ecutwfc"))
    rho0 = qein.get(parsed, "system", "ecutrho")
    ratio = float(rho0) / ecut0 if rho0 else 4.0
    kinfo = qein.kpoints_info(parsed)
    atoms = qein.read_structure(text)
    if kind == "ecut":
        vals = [float(v) for v in re.split(r"[,\s]+", str(opts.get("ecut_list") or "")) if v] or \
            [round(ecut0 * f) for f in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0)]
    else:
        if not kinfo.get("grid"):
            job.step("conv_k", status="skipped", error="K_POINTS is not 'automatic'")
            return
        vals = [int(v) for v in re.split(r"[,\s]+", str(opts.get("k_list") or "")) if v] or \
            sorted({max(1, int(round(max(kinfo["grid"]) * f))) for f in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5)})
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
            grid = qein.kgrid_for_density(atoms, v)
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
    gaps = [summary.get(k) for k in ("gap_bands", "gap_nscf", "gap_scf_grid") if summary.get(k)]
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
    # checks
    checks = []
    for out in sorted(d.glob("*.out")) + sorted(d.glob("conv_*/*/scf.out")):
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
    if res.get("projwfc", {}).get("spilling") and res["projwfc"]["spilling"] > 0.05:
        checks.append({"level": "info", "msg": f"PDOS spilling parameter {res['projwfc']['spilling']:.3f}: projections are approximate."})
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
                s = json.loads(r.read_text()).get("summary", {})
                item.update(energy_ry=s.get("energy_ry"), gap_ev=s.get("gap_ev"), is_metal=s.get("is_metal"))
            except json.JSONDecodeError:
                pass
        out.append(item)
    return out
