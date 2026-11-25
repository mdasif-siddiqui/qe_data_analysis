"""QE Studio web server. Run inside WSL:  bash start.sh  (then open http://localhost:8050)."""
from __future__ import annotations

import io
import json
import os
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

import numpy as np
from flask import Flask, abort, jsonify, request, send_file, send_from_directory

from qestudio import cif, config, jobs, nmr, phonon, pseudos, qein, qeout

ROOT = Path(__file__).resolve().parent
app = Flask(__name__, static_folder=str(ROOT / "static"), static_url_path="/static")
app.config["MAX_CONTENT_LENGTH"] = config.MAX_UPLOAD_MB * 1024 * 1024


@app.before_request
def _auth():
    """Optional shared password (HTTP basic auth; any user name)."""
    if not config.PASSWORD:
        return None
    import hmac
    a = request.authorization
    if a and a.password and hmac.compare_digest(a.password, config.PASSWORD):
        return None
    from flask import Response
    return Response("Password required", 401, {"WWW-Authenticate": 'Basic realm="QE Studio"'})


def local_only():
    if config.PUBLIC:
        abort(403, "Not available on the public server.")


@app.errorhandler(jobs.LimitError)
@app.errorhandler(ValueError)
def _limit(e):
    return jsonify({"error": str(e)}), 400


@app.errorhandler(413)
def _too_big(e):
    return jsonify({"error": f"Upload too large (limit {config.MAX_UPLOAD_MB} MB)."}), 413
_STATUS_CACHE: dict = {}


def win_path(p: Path) -> str:
    try:
        return subprocess.run(["wslpath", "-w", str(p)], capture_output=True, text=True, timeout=5).stdout.strip() or str(p)
    except OSError:
        return str(p)


def linux_path(p: str) -> str:
    p = p.strip().strip('"')
    if re.match(r"^[A-Za-z]:[\\/]", p) or p.startswith("\\\\"):
        try:
            return subprocess.run(["wslpath", "-u", p], capture_output=True, text=True, timeout=5).stdout.strip()
        except OSError:
            pass
    return os.path.expanduser(p)


def job_or_404(jid: str) -> jobs.Job:
    if not re.fullmatch(r"[\w.-]+", jid) or not (jobs.RUNS / jid).is_dir():
        abort(404)
    return jobs.manager().get(jid)


# ----------------------------------------------------------------------------- pages
@app.get("/")
def index():
    return send_from_directory(ROOT / "static", "index.html")


@app.get("/vendor/plotly.min.js")
def plotly_js():
    import plotly
    return send_from_directory(Path(plotly.__file__).parent / "package_data", "plotly.min.js")


# ----------------------------------------------------------------------------- status
@app.get("/api/status")
def status():
    if "qe_version" not in _STATUS_CACHE:
        _STATUS_CACHE["qe_version"] = jobs.qe_version()
    m = jobs.manager()
    return jsonify({
        "qe_version": _STATUS_CACHE["qe_version"],
        "pw_path": jobs.exe("pw.x"),
        "logical_cpus": os.cpu_count(),
        "physical_cores": jobs.physical_cores(),
        "launcher": m.launcher,
        "running": m.current.id if m.current else None,
        "queued": m.q.qsize(),
        "runs_dir": "" if config.PUBLIC else win_path(jobs.RUNS),
        "pseudo_library": pseudos.library_listing(),
        "public": config.PUBLIC,
        "gipaw": nmr.available(),
        "limits": {"atoms": config.MAX_ATOMS, "ecutwfc": config.MAX_ECUTWFC, "kpoints": config.MAX_KPOINTS,
                   "step_seconds": config.MAX_STEP_SECONDS, "keep_days": config.KEEP_DAYS} if config.PUBLIC else None,
    })


EXAMPLES = {"si": "si.scf.in", "al": "al.scf.in"}


@app.get("/api/example")
def example():
    name = EXAMPLES.get(request.args.get("name", "si"), "si.scf.in")
    return jsonify({"name": name, "text": (ROOT / "examples" / name).read_text()})


# ----------------------------------------------------------------------------- inspect / submit
def _read_upload():
    text, name = None, "scf.in"
    f = request.files.get("scf")
    if f:
        text = f.read().decode("utf-8", "replace")
        name = f.filename or name
    elif request.form.get("text"):
        text = request.form["text"]
        name = request.form.get("name") or name
    if text is None:
        abort(400, "No input file given")
    pp = {}
    for pf in request.files.getlist("pseudos"):
        if pf and pf.filename:
            pp[os.path.basename(pf.filename)] = pf.read()
    return text, name, pp


@app.post("/api/inspect")
def inspect():
    text, name, pp = _read_upload()
    res = qein.inspect(text)
    res["name"] = name
    if "namelists" in res:
        try:
            jobs.check_limits(qein.normalize_text(text)[0], qein.parse_input(qein.normalize_text(text)[0]), {})
        except jobs.LimitError as e:
            res.setdefault("errors", []).append(str(e))
    # pseudopotential availability
    avail = []
    for s in res.get("species", []):
        n = s["pseudo"]
        try:
            pseudos.safe_name(n)
        except ValueError as e:
            res.setdefault("errors", []).append(str(e))
            continue
        if n in pp:
            src = "uploaded"
        else:
            found = pseudos.find(n)
            if found and found.parent == pseudos.LIBRARY:
                src = "in QE Studio's pseudo library"
            elif found:
                src = f"found in {found.parent}"
            else:
                src = "will be downloaded from the QE library"
        avail.append({"name": n, "source": src, "label": s["label"]})
    res["pseudo_status"] = avail
    res["physical_cores"] = jobs.physical_cores()
    return jsonify(jobs._jsonable(res))


@app.post("/api/convert-cif")
def convert_cif():
    f = request.files.get("cif")
    text = f.read().decode("utf-8", "replace") if f else request.form.get("text", "")
    if not text.strip():
        abort(400, "No CIF given")
    try:
        r = cif.cif_to_scf(text, functional=request.form.get("functional", "pbe"),
                           primitive=request.form.get("primitive", "1") == "1",
                           occupations=request.form.get("occupations", "smearing"),
                           kspacing=float(request.form.get("kspacing") or 0.25))
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    return jsonify(jobs._jsonable(r))


@app.post("/api/jobs")
def submit():
    text, name, pp = _read_upload()
    check = qein.inspect(text)
    if check.get("errors"):
        return jsonify({"error": "Input has errors: " + " ".join(check["errors"])}), 400
    try:
        opts = json.loads(request.form.get("options") or "{}")
    except json.JSONDecodeError:
        opts = {}
    job = jobs.create_job(text, name, pp, opts)
    return jsonify({"id": job.id})


@app.get("/api/jobs")
def list_jobs():
    items = jobs.list_jobs()
    if config.PUBLIC:  # only the runs this browser created (ids sent by the page)
        mine = set(filter(None, request.args.get("ids", "").split(",")))
        items = [j for j in items if j["id"] in mine]
    return jsonify(jobs._jsonable(items))


@app.get("/api/jobs/<jid>")
def job_detail(jid):
    job = job_or_404(jid)
    st = job.load()
    log = job.dir / "pipeline.log"
    st["log"] = log.read_text()[-8000:] if log.exists() else ""
    # live tail of the currently running program
    cur = st.get("current_step")
    names = {"scf": "scf.out", "bands": "bands.out", "nscf": "nscf.out", "dos": "dos.out", "pdos": "projwfc.out",
             "ph_scf": "phonon/scf.out", "ph": "phonon/ph.out", "q2r": "phonon/q2r.out", "matdyn": "phonon/matdyn_dos.out",
             "sp_scf": "spectra/scf.out", "sp_ph": "spectra/ph.out", "dynmat": "spectra/dynmat.out",
             "uv_nscf": "uvvis/nscf.out", "nmr_scf": "nmr/scf.out", "nmr": "nmr/gipaw.out",
             "bandsx": "bands_pp.out"}
    if st.get("status") == "running" and cur in names and (job.dir / names[cur]).exists():
        st["live"] = "\n".join((job.dir / names[cur]).read_text(errors="replace").splitlines()[-25:])
    res = job.dir / "results.json"
    if res.exists() and st.get("status") not in ("running", "queued"):
        try:
            st["results"] = json.loads(res.read_text())
        except json.JSONDecodeError:
            pass
    st["files"] = sorted(str(p.relative_to(job.dir)) for p in job.dir.rglob("*") if p.is_file())
    st["win_dir"] = "" if config.PUBLIC else win_path(job.dir)
    st["public"] = config.PUBLIC
    return jsonify(jobs._jsonable(st))


@app.post("/api/jobs/<jid>/cancel")
def cancel(jid):
    job_or_404(jid)
    jobs.manager().cancel(jid)
    return jsonify({"ok": True})


@app.post("/api/jobs/<jid>/reanalyse")
def reanalyse(jid):
    job = job_or_404(jid)
    jobs.analyse_job(job)
    # after a successful re-analysis the run's status follows its calculation steps
    job.step("analysis", status="done")
    st = job.load()
    if st.get("status") not in ("running", "queued", "cancelled"):
        failed = [x["name"] for x in st["steps"] if x.get("status") == "failed"]
        job.update(status="done" if not failed else "done_with_errors", failed_steps=failed, error=None, current_step=None)
    return jsonify({"ok": True})


@app.delete("/api/jobs/<jid>")
def delete(jid):
    local_only()
    job = job_or_404(jid)
    if job.load().get("status") in ("running", "queued"):
        return jsonify({"error": "Cancel the job before deleting it."}), 400
    shutil.rmtree(job.dir)
    shutil.rmtree(job.scratch, ignore_errors=True)
    return jsonify({"ok": True})


@app.get("/api/jobs/<jid>/file/<path:rel>")
def job_file(jid, rel):
    job = job_or_404(jid)
    p = (job.dir / rel).resolve()
    if job.dir.resolve() not in p.parents or not p.is_file():
        abort(404)
    return send_file(p, as_attachment=request.args.get("dl") == "1", download_name=p.name,
                     mimetype="text/plain" if p.suffix in (".in", ".out", ".err", ".log", ".dat", ".gnu", ".dos") or "pdos" in p.name else None)


@app.get("/api/jobs/<jid>/zip")
def job_zip(jid):
    job = job_or_404(jid)
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for p in job.dir.rglob("*"):
            if p.is_file():
                z.write(p, f"{jid}/{p.relative_to(job.dir)}")
    buf.seek(0)
    return send_file(buf, as_attachment=True, download_name=f"{jid}.zip", mimetype="application/zip")


@app.post("/api/jobs/<jid>/open")
def open_folder(jid):
    local_only()
    job = job_or_404(jid)
    try:
        subprocess.Popen(["explorer.exe", win_path(job.dir)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except OSError as e:
        return jsonify({"error": str(e)}), 500
    return jsonify({"ok": True})


# ----------------------------------------------------------------------------- analyse existing outputs
def _phonon_file(p: Path) -> dict | None:
    """matdyn.x frequency (.freq) or phonon DOS files; None if p is neither."""
    try:
        head = p.read_text(errors="replace")[:300]
    except OSError:
        return None
    if "&plot" in head and "nbnd" in head:
        f = phonon.parse_freq(p)
        x = qeout.kpath_distance(f["q"])
        ticks = []
        for i in qeout.detect_kinks(f["q"]):
            lbl = "Γ" if np.allclose(f["q"][i], 0, atol=1e-6) else ""
            ticks.append({"x": float(x[i]), "label": lbl})
        fmin = float(f["freq"].min())
        return {"file": p.name, "program": "MATDYN", "calculation": "phonon dispersion", "job_done": None, "errors": [],
                "phonon_freq": {"x": x.tolist(), "freq": f["freq"].T.tolist(), "ticks": ticks,
                                "min_cm1": fmin, "max_cm1": float(f["freq"].max()),
                                "stable": fmin > phonon.IMAG_THRESHOLD_CM1}}
    if re.search(r"Frequency\s*\[cm\^-1\]", head):
        d = phonon.parse_phdos(p)
        natoms = max(1, len(d["pdos_atoms"]))
        th = phonon.thermodynamics(d["freq"], d["dos"], natoms)
        return {"file": p.name, "program": "MATDYN", "calculation": "phonon DOS", "job_done": None, "errors": [],
                "phonon_dos": {"freq": d["freq"].tolist(), "total": d["dos"].tolist(), "natoms": natoms,
                               "natoms_from_pdos": len(d["pdos_atoms"]) > 0, "thermo": th,
                               "zpe_mev_atom": th["zpe_ev_atom"] * 1000, "cv_300": phonon.value_at(th, "Cv", 300)}}
    return None


def _ph_out_summary(d: dict, path: Path) -> None:
    """Add the lowest/highest frequency over all q-points of a ph.x output."""
    txt = path.read_text(errors="replace")
    fr = [float(x) for x in re.findall(r"freq \(\s*\d+\)\s*=\s*[-\d.]+\s*\[THz\]\s*=\s*([-\d.]+)\s*\[cm-1\]", txt)]
    if fr:
        d["ph_min_cm1"], d["ph_max_cm1"] = min(fr), max(fr)
        d["ph_nq"] = len(re.findall(r"Diagonalizing the dynamical matrix", txt))


def _analyse_paths(paths: list[Path]) -> dict:
    outs = [p for p in paths if p.suffix.lower() in (".out", ".log", ".txt") or ".out" in p.name]
    ins = {p.name: p for p in paths if p.suffix.lower() == ".in"}
    items = []
    for p in sorted(paths):
        if p in outs or p.suffix.lower() in (".in", ".upf", ".gp", ".json", ".err"):
            continue
        try:
            ph = _phonon_file(p)
        except Exception as e:
            ph = {"file": p.name, "program": "MATDYN", "job_done": None, "errors": [f"could not parse: {e}"]}
        if ph:
            items.append(ph)
    for p in sorted(outs):
        try:
            d = qeout.parse_any(p)
        except Exception as e:
            items.append({"file": p.name, "program": "?", "errors": [f"could not parse: {e}"]})
            continue
        if d.get("program") == "unknown":
            continue
        if d.get("program") == "PHONON":
            _ph_out_summary(d, p)
        items.append(d)
    # reference levels (VBM / E_F) from SCF files for band plots
    scfs = [d for d in items if d.get("program") == "PWSCF" and d.get("calculation") in ("scf", "relax", "vc-relax")]
    for d in items:
        if d.get("program") != "PWSCF":
            continue
        bd = qeout.bands_from_pw(d)
        d.pop("bands", None)
        if not bd or len(bd["kpts"]) < 2:
            continue
        ref_src = next((s for s in scfs if s.get("nelec") == d.get("nelec") and s.get("nat") == d.get("nat")), None)
        fermi = d.get("fermi_ev") or (ref_src or {}).get("fermi_ev")
        fixed = fermi is None or d.get("homo_ev") is not None or (ref_src or {}).get("homo_ev") is not None
        edges = qeout.band_edges(bd["eig"], bd["kpts"], d.get("nelec"), fermi, fixed)
        d["gap"] = jobs._jsonable(edges)
        if d.get("calculation") == "bands":
            labels = _labels_from_inputs(ins, d, len(bd["kpts"]))
            if labels:
                idx, names, jumps = labels
                x = qeout.kpath_distance(bd["kpts"], jumps)
                ticks = []
                for i, n in zip(idx, names):
                    if i < len(x):
                        if ticks and abs(ticks[-1]["x"] - x[i]) < 1e-9:
                            ticks[-1]["label"] += "|" + n
                        else:
                            ticks.append({"x": float(x[i]), "label": n})
            else:
                x = qeout.kpath_distance(bd["kpts"])
                ticks = [{"x": float(x[i]), "label": ""} for i in qeout.detect_kinks(bd["kpts"])]
            if edges.get("metal") is False:
                for w in ("vbm", "cbm"):
                    edges[f"{w}_where"] = jobs._where(x[edges[f"{w}_index"]], ticks)
            ref = edges.get("vbm_ev") if edges.get("metal") is False else fermi
            d["band_plot"] = {"x": x.tolist(), "eig": [e.T.tolist() for e in bd["eig"]], "ticks": ticks,
                              "reference_ev": ref, "reference_label": "VBM" if edges.get("metal") is False else "E_F",
                              "edges": jobs._jsonable(edges), "reference_from": (ref_src or {}).get("file")}
    return {"items": jobs._jsonable(items)}


def _labels_from_inputs(ins: dict, d: dict, nks: int):
    """Find a bands .in with a crystal_b path giving exactly nks points; read labels from comments."""
    cands = sorted(ins.values(), key=lambda p: (Path(d["file"]).stem not in p.stem, p.name))
    for p in cands:
        raw = p.read_text(errors="replace").replace("\r", "")
        m = re.search(r"(?im)^\s*K_POINTS\s*[{(]?\s*(crystal_b|tpiba_b)\s*[})]?\s*\n\s*(\d+)\s*\n", raw)
        if not m:
            continue
        n = int(m.group(2))
        lines = raw[m.end():].split("\n")[:n]
        w, names = [], []
        for ln in lines:
            parts = ln.split("!")
            nums = parts[0].split()
            if len(nums) < 4:
                break
            w.append(int(float(nums[3])))
            lbl = parts[1].strip() if len(parts) > 1 else (nums[4] if len(nums) > 4 else "")
            names.append(qein.pretty_label(lbl.strip()) if lbl else "")
        if len(w) != n or sum(w[:-1]) + 1 != nks:
            continue
        idx = [sum(w[:i]) for i in range(n)]
        jumps = [idx[i] for i in range(n - 1) if w[i] == 1]
        return idx, names, jumps
    return None


@app.post("/api/analyze")
def analyze():
    folder = request.form.get("folder", "").strip()
    if folder:
        local_only()
        p = Path(linux_path(folder))
        if not p.is_dir():
            return jsonify({"error": f"Folder not found: {folder} ({p})"}), 400
        paths = [q for q in p.iterdir() if q.is_file()]
        res = _analyse_paths(paths)
        res["folder"] = str(p)
        return jsonify(res)
    files = request.files.getlist("files")
    if not files:
        return jsonify({"error": "No files uploaded"}), 400
    with tempfile.TemporaryDirectory() as td:
        paths = []
        for f in files:
            name = os.path.basename(f.filename or "file.out")
            q = Path(td) / name
            f.save(q)
            paths.append(q)
        res = _analyse_paths(paths)
    for it in res["items"]:
        it.pop("path", None)
    return jsonify(res)


if __name__ == "__main__":
    jobs.manager()
    print(f"QE Studio running at http://localhost:{config.PORT}" + (" (public mode)" if config.PUBLIC else ""))
    app.run(host=config.HOST, port=config.PORT, debug=False, threaded=True)
