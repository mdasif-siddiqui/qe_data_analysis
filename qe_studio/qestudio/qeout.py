"""Parsers for pw.x, dos.x, projwfc.x and bands.x outputs, plus band-gap analysis."""
from __future__ import annotations

import glob
import math
import os
import re
from pathlib import Path

import numpy as np

RY_EV = 13.605693122994
FLOAT = r"[-+]?\d*\.\d+(?:[eEdD][-+]?\d+)?|[-+]?\d+(?:[eEdD][-+]?\d+)?"


def _f(s: str) -> float:
    return float(s.lower().replace("d", "e"))


def _last(pattern: str, text: str, group: int = 1, conv=_f):
    m = None
    for m in re.finditer(pattern, text, re.M):
        pass
    if m is None:
        return None
    try:
        return conv(m.group(group))
    except (ValueError, TypeError):
        return None


def _wall_seconds(s: str) -> float | None:
    s = s.strip()
    total, found = 0.0, False
    for val, unit in re.findall(r"([\d.]+)\s*([hms])", s):
        found = True
        total += float(val) * {"h": 3600, "m": 60, "s": 1}[unit]
    return total if found else None


def program_of(text: str) -> str:
    m = re.search(r"Program\s+(\S+)\s+v\.?\s*([\w.]+)", text)
    return m.group(1) if m else "unknown"


def errors_of(text: str) -> list[str]:
    errs = []
    for m in re.finditer(r"%{20,}\s*\n(.*?)\n\s*%{20,}", text, re.S):
        block = " ".join(l.strip() for l in m.group(1).splitlines() if l.strip())
        errs.append(block)
    return errs


# ----------------------------------------------------------------------------- pw.x
def _parse_eigen_block(block: str) -> dict:
    """Parse the k-point/eigenvalue listing that follows 'End of ... calculation'."""
    spins = []
    parts = re.split(r"-{3,}\s*SPIN\s+(UP|DOWN)\s*-{3,}", block)
    if len(parts) > 1:
        chunks = [parts[i + 1] for i in range(1, len(parts) - 1, 2)]
    else:
        chunks = [block]
    for chunk in chunks:
        kpts, eigs = [], []
        cur = None
        mode = None
        for line in chunk.splitlines():
            if "bands (ev)" in line or "band energies (ev)" in line:
                head = line.split("(")[0] if "PWs" in line else line.split("band")[0]
                nums = re.findall(r"-?\d+\.\d+", head.split("=", 1)[1]) if "=" in head else []
                kpts.append([float(x) for x in nums[:3]])
                cur = []
                eigs.append(cur)
                mode = "e"
                continue
            if "occupation numbers" in line:
                mode = None
                continue
            if re.search(r"Fermi|highest|lowest|total energy|Writing|SPIN", line):
                mode = None
                if cur is not None and re.search(r"Fermi|highest|total energy|Writing", line):
                    break
                continue
            if mode == "e":
                vals = re.findall(r"-?\d+\.\d+", line)
                if vals:
                    cur.extend(float(v) for v in vals)
        if eigs:
            nb = min(len(e) for e in eigs)
            spins.append({"kpts": kpts, "eig": [e[:nb] for e in eigs]})
    return {"spins": spins}


def parse_pw(text: str) -> dict:
    d: dict = {"program": "PWSCF"}
    m = re.search(r"Program PWSCF v\.?\s*([\w.]+)", text)
    d["version"] = m.group(1) if m else None
    d["job_done"] = "JOB DONE" in text
    d["errors"] = errors_of(text)
    d["not_converged"] = "convergence NOT achieved" in text
    # type of calculation
    if "new unit-cell volume" in text or "Variable-cell" in text:
        calc = "vc-relax"
    elif "BFGS Geometry Optimization" in text:
        calc = "relax"
    elif "Molecular Dynamics Calculation" in text:
        calc = "md"
    elif "End of band structure calculation" in text or "Band Structure Calculation" in text:
        calc = "nscf" if re.search(r"the Fermi energy|highest occupied", text) else "bands"
    else:
        calc = "scf"
    d["calculation"] = calc
    d["alat_bohr"] = _last(r"lattice parameter \(alat\)\s*=\s*(" + FLOAT + ")", text)
    d["volume_bohr3"] = _last(r"unit-cell volume\s*=\s*(" + FLOAT + ")", text)
    d["nat"] = _last(r"number of atoms/cell\s*=\s*(\d+)", text, conv=int)
    d["ntyp"] = _last(r"number of atomic types\s*=\s*(\d+)", text, conv=int)
    d["nelec"] = _last(r"number of electrons\s*=\s*(" + FLOAT + ")", text)
    d["nbnd"] = _last(r"number of Kohn-Sham states\s*=\s*(\d+)", text, conv=int)
    d["ecutwfc"] = _last(r"kinetic-energy cutoff\s*=\s*(" + FLOAT + ")", text)
    d["ecutrho"] = _last(r"charge density cutoff\s*=\s*(" + FLOAT + ")", text)
    d["nks"] = _last(r"number of k points\s*=\s*(\d+)", text, conv=int)
    d["nspin"] = 2 if re.search(r"SPIN UP|spin up/dw Fermi|LSDA", text) else 1
    d["mpi_procs"] = _last(r"Number of MPI processes:\s*(\d+)", text, conv=int) or \
        _last(r"running on\s+(\d+)\s+processor", text, conv=int)
    d["xc"] = _last(r"Exchange-correlation\s*=\s*(.+?)\s*$", text, conv=str)
    pp = re.findall(r"PseudoPot\.\s*#\s*\d+ for (\w+)\s+read from file:\s*\n\s*(\S+)", text)
    d["pseudos"] = [{"element": e, "file": os.path.basename(f)} for e, f in pp]

    # energies
    d["energy_ry"] = _last(r"^!\s+total energy\s*=\s*(" + FLOAT + ")", text)
    final = _last(r"Final (?:energy|enthalpy)\s*=\s*(" + FLOAT + ")", text)
    if final is not None:
        d["final_energy_ry"] = final
    for key, pat in (("one_electron", "one-electron contribution"), ("hartree", "hartree contribution"),
                     ("xc_energy", "xc contribution"), ("ewald", "ewald contribution"),
                     ("smearing_contrib", r"smearing contrib\. \(-TS\)")):
        v = _last(pat + r"\s*=\s*(" + FLOAT + ")", text)
        if v is not None:
            d.setdefault("energy_terms_ry", {})[key] = v
    d["fermi_ev"] = _last(r"the Fermi energy is\s+(" + FLOAT + ")", text)
    m = None
    for m in re.finditer(r"the spin up/dw Fermi energies are\s+(" + FLOAT + r")\s+(" + FLOAT + ")", text):
        pass
    if m:
        d["fermi_up_dw_ev"] = [_f(m.group(1)), _f(m.group(2))]
    m = None
    for m in re.finditer(r"highest occupied, lowest unoccupied level \(ev\):\s+(" + FLOAT + r")\s+(" + FLOAT + ")", text):
        pass
    if m:
        d["homo_ev"], d["lumo_ev"] = _f(m.group(1)), _f(m.group(2))
        d["homo_lumo_gap_ev"] = d["lumo_ev"] - d["homo_ev"]
    else:
        d["homo_ev"] = _last(r"highest occupied level \(ev\):\s+(" + FLOAT + ")", text)
    d["total_force"] = _last(r"Total force\s*=\s*(" + FLOAT + ")", text)
    d["pressure_kbar"] = _last(r"P=\s*(" + FLOAT + ")", text)
    d["total_mag"] = _last(r"total magnetization\s*=\s*(" + FLOAT + ")", text)
    d["abs_mag"] = _last(r"absolute magnetization\s*=\s*(" + FLOAT + ")", text)
    # moments integrated in a small sphere per atom (last block printed); antiferromagnets have total ~ 0
    blocks = text.split("Magnetic moment per site")
    if len(blocks) > 1:
        sites = re.findall(r"atom\s*:?\s*(\d+)\s*(?:\(R=\s*(?:" + FLOAT + r")\))?\s*charge\s*[=:]\s*(" + FLOAT + r")\s*magn\s*[=:]\s*("
                           + FLOAT + ")", blocks[-1][:4000])
        if sites:
            d["site_moments"] = [{"atom": int(a), "charge": _f(c), "magn": _f(m)} for a, c, m in sites]
    d["scf_iterations"] = _last(r"convergence has been achieved in\s+(\d+)\s+iterations", text, conv=int)
    wall = _last(r"PWSCF\s*:\s*(.*?)CPU\s+(.*?)WALL", text, group=2, conv=str)
    d["wall_s"] = _wall_seconds(wall) if wall else None

    # SCF history, split into ionic steps at every '!' line
    cycles, cur = [], []
    for line in text.splitlines():
        mi = re.match(r"\s*iteration #\s*(\d+)", line)
        if mi:
            cur.append({"iter": int(mi.group(1)), "energy": None, "accuracy": None})
            continue
        me = re.match(r"\s+total energy\s*=\s*(" + FLOAT + ")", line)
        if me and cur:
            cur[-1]["energy"] = _f(me.group(1))
            continue
        ma = re.match(r"\s*estimated scf accuracy\s*<\s*(" + FLOAT + ")", line)
        if ma and cur:
            cur[-1]["accuracy"] = _f(ma.group(1))
            continue
        if line.startswith("!") and cur:
            mm = re.search(r"=\s*(" + FLOAT + ")", line)
            if mm:
                cur[-1]["energy"] = _f(mm.group(1))
            cycles.append(cur)
            cur = []
    if cur:
        cycles.append(cur)
    # the converged iteration prints its accuracy after '!'; recover it
    accs = [_f(x) for x in re.findall(r"estimated scf accuracy\s*<\s*(" + FLOAT + ")", text)]
    if cycles and accs and cycles[-1] and cycles[-1][-1]["accuracy"] is None:
        cycles[-1][-1]["accuracy"] = accs[-1]
    d["scf_cycles"] = cycles

    # relaxation trajectory
    if calc in ("relax", "vc-relax", "md"):
        energies = [_f(x) for x in re.findall(r"^!\s+total energy\s*=\s*(" + FLOAT + ")", text, re.M)]
        forces = [_f(x) for x in re.findall(r"Total force\s*=\s*(" + FLOAT + ")", text)]
        press = [_f(x) for x in re.findall(r"P=\s*(" + FLOAT + ")", text)]
        vols = [_f(x) for x in re.findall(r"new unit-cell volume\s*=\s*(" + FLOAT + ")", text)]
        d["relax"] = {"energy_ry": energies, "force": forces, "pressure": press, "volume_bohr3": vols,
                      "converged": bool(re.search(r"bfgs converged|End of BFGS Geometry Optimization", text))
                      and "The maximum number of steps has been reached" not in text}

    # eigenvalues
    idx = max(text.rfind("End of self-consistent calculation"), text.rfind("End of band structure calculation"))
    if idx >= 0:
        bands = _parse_eigen_block(text[idx:])
        if bands["spins"]:
            d["bands"] = bands
    return d


# ----------------------------------------------------------------------------- band analysis
def kpath_distance(kpts: np.ndarray, jumps_after: list[int] | None = None) -> np.ndarray:
    """Cumulative distance along k (units of the input, i.e. 2pi/alat); zero step after each jump index."""
    dk = np.linalg.norm(np.diff(kpts, axis=0), axis=1)
    for j in jumps_after or []:
        if 0 <= j < len(dk):
            dk[j] = 0.0
    return np.concatenate([[0.0], np.cumsum(dk)])


def detect_kinks(kpts: np.ndarray, tol: float = 1e-3) -> list[int]:
    """Indices where the path changes direction (or jumps) – used when labels are unknown."""
    if len(kpts) < 3:
        return [0, len(kpts) - 1]
    idx = [0]
    d = np.diff(kpts, axis=0)
    n = np.linalg.norm(d, axis=1)
    typical = np.median(n[n > 1e-8]) if np.any(n > 1e-8) else 0
    for i in range(1, len(d)):
        a, b = d[i - 1], d[i]
        if n[i - 1] < 1e-8 or n[i] < 1e-8:
            continue
        cos = a @ b / (n[i - 1] * n[i])
        if cos < 1 - tol or (typical and n[i] > 4 * typical):
            idx.append(i)
    idx.append(len(kpts) - 1)
    return sorted(set(idx))


def band_edges(eigs: list[np.ndarray], kpts: np.ndarray, nelec: float | None, fermi: float | None,
               occupations_fixed: bool, noncolin: bool = False) -> dict:
    """VBM/CBM/gap from eigenvalues of one or two spin channels.

    eigs: list (per spin) of arrays (nk, nbnd) in eV.
    """
    res: dict = {}
    all_e = np.concatenate([e.ravel() for e in eigs])
    nspin = len(eigs)
    vb_list, cb_list = [], []
    if nelec is not None and nspin == 1 and (occupations_fixed or fermi is None):
        nocc = nelec if noncolin else nelec / 2
        if abs(nocc - round(nocc)) < 1e-6 and 0 < round(nocc) < eigs[0].shape[1]:
            n = int(round(nocc))
            vb = eigs[0][:, n - 1]
            cb = eigs[0][:, n]
            res["method"] = f"band count (VB = band {n}, CB = band {n+1})"
            vb_list, cb_list = [(vb, 0)], [(cb, 0)]
    if not vb_list:
        ref = fermi
        if ref is None:
            return {"metal": None, "method": "no Fermi level or electron count available"}
        res["method"] = "relative to the Fermi level"
        for s, e in enumerate(eigs):
            below = np.where(e <= ref, e, -np.inf)
            above = np.where(e > ref, e, np.inf)
            # a band crossing E_F means metal
            for b in range(e.shape[1]):
                if e[:, b].min() < ref - 1e-4 and e[:, b].max() > ref + 1e-4:
                    res.update(metal=True, gap_ev=0.0, crossing_band=b + 1)
                    res["vbm_ev"] = ref
                    res["cbm_ev"] = ref
                    res["reference_ev"] = ref
                    return res
            vb_list.append((below.max(axis=1), s))
            cb_list.append((above.min(axis=1), s))
    vbm, vbm_k, vbm_s = -np.inf, 0, 0
    for arr, s in vb_list:
        i = int(np.argmax(arr))
        if arr[i] > vbm:
            vbm, vbm_k, vbm_s = float(arr[i]), i, s
    cbm, cbm_k, cbm_s = np.inf, 0, 0
    for arr, s in cb_list:
        i = int(np.argmin(arr))
        if arr[i] < cbm:
            cbm, cbm_k, cbm_s = float(arr[i]), i, s
    if not np.isfinite(vbm) or not np.isfinite(cbm):
        return {"metal": None, "method": "not enough bands to find the conduction band"}
    gap = cbm - vbm
    if gap <= 1e-4:
        return {"metal": True, "gap_ev": 0.0, "vbm_ev": vbm, "cbm_ev": cbm, "reference_ev": fermi if fermi is not None else vbm,
                "method": res.get("method")}
    direct = np.inf
    for (vb, s1), (cb, s2) in zip(vb_list, cb_list):
        direct = min(direct, float(np.min(cb - vb)))
    res.update(
        metal=False, gap_ev=gap, vbm_ev=vbm, cbm_ev=cbm, reference_ev=vbm,
        vbm_k=kpts[vbm_k].tolist(), cbm_k=kpts[cbm_k].tolist(), vbm_index=vbm_k, cbm_index=cbm_k,
        vbm_spin=vbm_s, cbm_spin=cbm_s,
        direct_gap_ev=direct,
        gap_type="direct" if np.allclose(kpts[vbm_k], kpts[cbm_k], atol=1e-4) else "indirect",
    )
    return res


# ----------------------------------------------------------------------------- dos.x / projwfc.x
def parse_dos_file(path: str | Path) -> dict:
    path = Path(path)
    head = path.open().readline()
    ef = re.search(r"EFermi\s*=\s*(" + FLOAT + ")", head)
    data = np.loadtxt(path, comments="#")
    out = {"fermi_ev": _f(ef.group(1)) if ef else None, "energy": data[:, 0].tolist()}
    if "dosup" in head:
        out["dos_up"] = data[:, 1].tolist()
        out["dos_dw"] = data[:, 2].tolist()
        out["idos"] = data[:, 3].tolist()
    else:
        out["dos"] = data[:, 1].tolist()
        out["idos"] = data[:, 2].tolist()
    return out


def parse_pdos(folder: str | Path, filpdos: str) -> dict:
    """Aggregate projwfc.x files per element and orbital (s, p, d, f)."""
    folder = Path(folder)
    groups: dict = {}
    energy = None
    spin = False
    for f in sorted(glob.glob(str(folder / f"{filpdos}.pdos_atm#*"))):
        m = re.search(r"pdos_atm#(\d+)\((\w+)\)_wfc#(\d+)\((\w)", os.path.basename(f))
        if not m:
            continue
        elem, orb = m.group(2), m.group(4)
        with open(f) as fh:
            head = fh.readline()
        data = np.loadtxt(f, comments="#")
        if data.ndim == 1:
            continue
        energy = data[:, 0]
        key = f"{elem} {orb}"
        if "ldosup" in head:
            spin = True
            up, dw = data[:, 1], data[:, 2]
            g = groups.setdefault(key, {"up": np.zeros_like(up), "dw": np.zeros_like(dw)})
            g["up"] += up
            g["dw"] += dw
        else:
            v = data[:, 1]
            g = groups.setdefault(key, {"total": np.zeros_like(v)})
            g["total"] += v
    out = {"energy": energy.tolist() if energy is not None else [], "spin": spin, "projections": {}}
    for k, g in sorted(groups.items()):
        out["projections"][k] = {kk: vv.tolist() for kk, vv in g.items()}
    tot = folder / f"{filpdos}.pdos_tot"
    if tot.exists():
        data = np.loadtxt(tot, comments="#")
        out["total_dos"] = data[:, 1].tolist() if not spin else (data[:, 1]).tolist()
        if spin:
            out["total_dos_dw"] = data[:, 2].tolist()
    return out


def parse_projwfc_out(text: str) -> dict:
    res = {"job_done": "JOB DONE" in text, "errors": errors_of(text), "lowdin": []}
    blk = text.split("Lowdin Charges:")[-1] if "Lowdin Charges:" in text else ""
    atoms: dict = {}
    for m in re.finditer(r"Atom #\s*(\d+): total charge =\s*(" + FLOAT + r")(.*)", blk):
        orb = dict((k, _f(v)) for k, v in re.findall(r"(\w+)\s*=\s*(" + FLOAT + ")", m.group(3)))
        a = atoms.setdefault(int(m.group(1)), {"atom": int(m.group(1)), "charge": _f(m.group(2))})
        a.update(orb)
    res["lowdin"] = [atoms[k] for k in sorted(atoms)]
    sp = re.search(r"Spilling Parameter:\s*(" + FLOAT + ")", text)
    res["spilling"] = _f(sp.group(1)) if sp else None
    return res


# ----------------------------------------------------------------------------- generic
def parse_any(path: str | Path) -> dict:
    """Summary of any QE output file (used by the 'Analyse outputs' page)."""
    path = Path(path)
    text = path.read_text(errors="replace")
    prog = program_of(text)
    base = {"file": path.name, "path": str(path), "program": prog,
            "job_done": "JOB DONE" in text, "errors": errors_of(text)}
    if prog == "PWSCF":
        d = parse_pw(text)
        d.update({k: v for k, v in base.items() if k not in d or k == "file"})
        d["file"] = path.name
        d["path"] = str(path)
        return d
    if prog == "PROJWFC":
        base.update(parse_projwfc_out(text))
    m = re.search(r"Program\s+\S+\s+v\.?\s*([\w.]+)", text)
    base["version"] = m.group(1) if m else None
    return base


def bands_from_pw(d: dict) -> dict | None:
    """Band data (eV) from a parsed pw.x output: {'kpts', 'eig': [per spin (nk,nb)]}."""
    b = d.get("bands")
    if not b or not b["spins"]:
        return None
    kpts = np.array(b["spins"][0]["kpts"], dtype=float)
    eig = [np.array(s["eig"], dtype=float) for s in b["spins"]]
    return {"kpts": kpts, "eig": eig}
