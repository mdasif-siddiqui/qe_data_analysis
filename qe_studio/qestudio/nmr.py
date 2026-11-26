"""NMR shielding with QE-GIPAW (gipaw.x), built from source by build_gipaw.sh into ~/gipaw."""
from __future__ import annotations

import re
from pathlib import Path

GIPAW_HOME = Path.home() / "gipaw"
GIPAW_BIN = GIPAW_HOME / "bin"
BUILD_ENV = Path.home() / "miniforge3" / "envs" / "qebuild"
F = r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?"

# NMR-active isotope usually measured for each element (for labels only)
ISOTOPES = {"H": "¹H", "Li": "⁷Li", "B": "¹¹B", "C": "¹³C", "N": "¹⁵N", "O": "¹⁷O", "F": "¹⁹F", "Na": "²³Na",
            "Mg": "²⁵Mg", "Al": "²⁷Al", "Si": "²⁹Si", "P": "³¹P", "S": "³³S", "Cl": "³⁵Cl", "K": "³⁹K", "Ca": "⁴³Ca",
            "Ti": "⁴⁷Ti", "V": "⁵¹V", "Zn": "⁶⁷Zn", "Ga": "⁷¹Ga", "Ge": "⁷³Ge", "As": "⁷⁵As", "Se": "⁷⁷Se",
            "Sn": "¹¹⁹Sn", "Pb": "²⁰⁷Pb", "Cd": "¹¹³Cd", "Y": "⁸⁹Y", "Zr": "⁹¹Zr", "Nb": "⁹³Nb", "Ba": "¹³⁷Ba"}


def available() -> bool:
    return (GIPAW_BIN / "gipaw.x").exists() and (GIPAW_BIN / "pw.x").exists()


def env() -> dict:
    """Environment for the self-built programs (their own MPI and libraries)."""
    import os
    return {"PATH": f"{GIPAW_BIN}:{BUILD_ENV / 'bin'}:{os.environ.get('PATH', '')}",
            "LD_LIBRARY_PATH": f"{BUILD_ENV / 'lib'}:{os.environ.get('LD_LIBRARY_PATH', '')}"}


def gipaw_ready(path: Path) -> tuple[bool, str]:
    """GIPAW needs reconstruction data: PAW files with paw_as_gipaw, or PP_GIPAW orbitals with non-zero radii."""
    head = path.read_text(errors="replace")
    if not re.search(r'has_gipaw="(true|T|\.true\.)"', head, re.I) and "<PP_GIPAW" not in head:
        return False, "no GIPAW data"
    if re.search(r'paw_as_gipaw="(true|T|\.true\.)"', head, re.I):
        return True, "PAW (paw_as_gipaw)"
    radii = [float(x) for x in re.findall(r'<PP_GIPAW_ORBITAL\.\d+[^>]*?cutoff_radius="([-+\d.eE]+)"', head)]
    if radii and min(radii) > 0:
        return True, "GIPAW orbitals"
    return False, "GIPAW data present but with zero cutoff radii (PSlibrary ultrasoft files) - use the PAW version"


def gipaw_input(prefix: str, tmp_dir: str, diagonalization: str | None = None) -> str:
    diag = f"  diagonalization = '{diagonalization}'\n" if diagonalization else ""
    return (f"&INPUTGIPAW\n  job = 'nmr'\n  prefix = '{prefix}'\n  tmp_dir = '{tmp_dir}'\n  q_gipaw = 0.01\n"
            f"  spline_ps = .true.\n  use_nmr_macroscopic_shape = .true.\n{diag}/\n")


def parse_gipaw(text: str) -> dict:
    """Shielding tensors per atom from the 'Total NMR chemical shifts' block of gipaw.x."""
    res: dict = {"job_done": "JOB DONE" in text, "atoms": []}
    blk = text.split("Total NMR chemical shifts in ppm", 1)
    if len(blk) < 2:
        return res
    b = blk[1]
    for m in re.finditer(r"Atom\s+(\d+)\s+(\w+)\s+pos:\s*\(\s*(" + F + r")\s+(" + F + r")\s+(" + F + r")\s*\)\s+Total sigma:\s*(" + F + r")\s*\n"
                         r"((?:\s*" + F + r"\s+" + F + r"\s+" + F + r"\s*\n){3})", b):
        idx, el = int(m.group(1)), m.group(2)
        tensor = [[float(x) for x in re.findall(F, l)] for l in m.group(7).strip().splitlines()]
        a = {"atom": idx, "element": el, "pos_alat": [float(m.group(i)) for i in (3, 4, 5)],
             "sigma_iso": float(m.group(6)), "tensor": tensor}
        an = re.search(rf"(?:^|\n)\s*{re.escape(el)}\s+{idx}\s+anisotropy:\s*({F})\s+eta:\s*({F})", b[m.end():m.end() + 600])
        if an:
            a["anisotropy"], a["eta"] = float(an.group(1)), float(an.group(2))
        pcs = re.findall(rf"\n\s*{re.escape(el)}\s+{idx}\s+sigma_(xx|yy|zz)=\s*({F})", b[m.end():m.end() + 900])
        if pcs:
            a["principal"] = {k: float(v) for k, v in pcs}
        res["atoms"].append(a)
    by_el: dict = {}
    for a in res["atoms"]:
        by_el.setdefault(a["element"], []).append(a["sigma_iso"])
    res["mean_sigma"] = {k: sum(v) / len(v) for k, v in by_el.items()}
    res["isotopes"] = {k: ISOTOPES.get(k, k) for k in by_el}
    m = re.search(r"f-sum rule \(should be\s*(" + F + r")\):\s*\n((?:\s*" + F + r"\s+" + F + r"\s+" + F + r"\s*\n){3})", text)
    if m:
        rows = [[float(x) for x in re.findall(F, l)] for l in m.group(2).strip().splitlines()]
        res["fsum_target"] = float(m.group(1))
        res["fsum_trace"] = sum(rows[i][i] for i in range(3)) / 3
    return res
