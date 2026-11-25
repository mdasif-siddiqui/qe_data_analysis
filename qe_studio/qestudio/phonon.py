"""Phonons: input generation for ph.x / q2r.x / matdyn.x, parsers, and harmonic thermodynamics."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

CM1_EV = 1.239841984e-4        # 1 cm^-1 in eV
CM1_THZ = 0.0299792458         # 1 cm^-1 in THz
KB_EV = 8.617333262e-5         # Boltzmann constant, eV/K
EV_JMOL = 96485.33212          # 1 eV per particle in J/mol
IMAG_THRESHOLD_CM1 = -10.0     # below this a frequency counts as imaginary (numerical noise is smaller)


# ----------------------------------------------------------------------------- inputs
def ph_input(prefix: str, outdir: str, q: list[int], tr2_ph: float = 1e-14) -> str:
    return (f"phonons\n&INPUTPH\n  prefix = '{prefix}'\n  outdir = '{outdir}'\n  fildyn = '{prefix}.dyn'\n"
            f"  tr2_ph = {tr2_ph:.1e}\n  ldisp = .true.\n  nq1 = {q[0]}, nq2 = {q[1]}, nq3 = {q[2]}\n/\n")


def q2r_input(prefix: str) -> str:
    return f"&INPUT\n  fildyn = '{prefix}.dyn'\n  zasr = 'crystal'\n  flfrc = '{prefix}.fc'\n/\n"


def matdyn_path_input(prefix: str, kpath_card: str) -> str:
    """Dispersion along the same crystal_b path as the band structure."""
    pts = []
    for ln in kpath_card.strip().splitlines()[2:]:
        nums = ln.split("!")[0].split()
        if len(nums) >= 4:
            pts.append(f"  {nums[0]} {nums[1]} {nums[2]} {nums[3]}")
    return (f"&INPUT\n  asr = 'crystal'\n  flfrc = '{prefix}.fc'\n  flfrq = '{prefix}.freq'\n"
            f"  q_in_band_form = .true.\n  q_in_cryst_coord = .true.\n/\n{len(pts)}\n" + "\n".join(pts) + "\n")


def matdyn_dos_input(prefix: str, grid: list[int], delta_cm: float = 1.0, degauss_cm: float = 3.0) -> str:
    """Phonon DOS with Gaussian broadening: matdyn's tetrahedron method can return NaN for every point."""
    return (f"&INPUT\n  asr = 'crystal'\n  flfrc = '{prefix}.fc'\n  flfrq = ''\n  dos = .true.\n"
            f"  fldos = '{prefix}.phdos'\n  nk1 = {grid[0]}, nk2 = {grid[1]}, nk3 = {grid[2]}\n"
            f"  deltaE = {delta_cm}\n  degauss = {degauss_cm}\n/\n")


def matdyn_mesh_input(prefix: str, grid: list[int]) -> str:
    """Frequencies on a full uniform q-mesh (crystal coordinates) for thermodynamics by direct summation."""
    pts = [f"  {i / grid[0]:.6f} {j / grid[1]:.6f} {k / grid[2]:.6f}"
           for i in range(grid[0]) for j in range(grid[1]) for k in range(grid[2])]
    return (f"&INPUT\n  asr = 'crystal'\n  flfrc = '{prefix}.fc'\n  flfrq = '{prefix}.mesh.freq'\n"
            f"  q_in_cryst_coord = .true.\n/\n{len(pts)}\n" + "\n".join(pts) + "\n")


# ----------------------------------------------------------------------------- parsers
def parse_freq(path: str | Path) -> dict:
    """matdyn.x flfrq file: q-points (cartesian, 2pi/alat) and frequencies (cm^-1)."""
    text = Path(path).read_text()
    m = re.search(r"nbnd\s*=\s*(\d+)\s*,\s*nks\s*=\s*(\d+)", text)
    if not m:
        raise ValueError(f"{path}: not a matdyn frequency file")
    nb, nq = int(m.group(1)), int(m.group(2))
    num = re.compile(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+")
    lines = [l for l in text[m.end():].splitlines() if num.search(l)]  # skips the closing ' /' of &plot
    qs, fr, i = [], [], 0
    # per q: one line with the coordinates (some QE versions add the path distance), then nb frequencies
    while i < len(lines) and len(qs) < nq:
        qs.append([float(x) for x in num.findall(lines[i])[:3]])
        i += 1
        vals: list[float] = []
        while len(vals) < nb and i < len(lines):
            vals += [float(x) for x in num.findall(lines[i])]
            i += 1
        fr.append(vals[:nb])
    if len(qs) != nq:
        raise ValueError(f"{path}: expected {nq} q-points, found {len(qs)}")
    return {"q": np.array(qs), "freq": np.array(fr)}


def parse_phdos(path: str | Path) -> dict:
    """matdyn.x DOS file: frequency (cm^-1), total DOS, and one projected DOS per atom."""
    data = np.loadtxt(path, comments="#")
    if np.isnan(data[:, 1]).all():
        raise ValueError(f"{Path(path).name}: the phonon DOS is NaN everywhere (matdyn tetrahedron method failed; use degauss)")
    data = np.nan_to_num(data)
    return {"freq": data[:, 0], "dos": data[:, 1], "pdos_atoms": data[:, 2:].T if data.shape[1] > 2 else np.empty((0, len(data)))}


def parse_dyn_gamma(path: str | Path) -> list[float]:
    """Raw Gamma frequencies (cm^-1) from fildyn1, before the acoustic sum rule is imposed."""
    text = Path(path).read_text()
    return [float(x) for x in re.findall(r"freq \(\s*\d+\)\s*=\s*[-\d.]+\s*\[THz\]\s*=\s*([-\d.]+)\s*\[cm-1\]", text)]


# ----------------------------------------------------------------------------- thermodynamics
def thermodynamics(freq_cm: np.ndarray, dos: np.ndarray, natoms: int,
                   temps: np.ndarray | None = None) -> dict:
    """Harmonic thermodynamics per mole of atoms from the phonon DOS.

    The DOS is normalised to 3*natoms modes. Imaginary (negative) frequencies are left out and
    their weight is reported, since they have no harmonic thermodynamics.
    """
    if temps is None:
        temps = np.arange(0.0, 1001.0, 10.0)
    w = np.asarray(freq_cm, float)
    g = np.clip(np.asarray(dos, float), 0, None)
    dw = np.gradient(w)
    total = np.sum(g * dw)
    neg = np.sum(g[w < 0] * dw[w < 0]) / total if total > 0 else 0.0
    pos = w > 1e-6
    w, g, dw = w[pos], g[pos], dw[pos]
    g = g * (3 * natoms) / max(total, 1e-30)          # normalise to 3N modes
    e = w * CM1_EV                                      # mode energies, eV
    zpe = np.sum(g * dw * e / 2)                         # eV per cell
    out = {"T": temps.tolist(), "Cv": [], "S": [], "F": [], "U": []}
    per_mol = EV_JMOL / natoms                           # eV per cell -> J per mole of atoms
    for T in temps:
        if T < 1e-8:
            out["Cv"].append(0.0)
            out["S"].append(0.0)
            out["F"].append(zpe * per_mol / 1000)
            out["U"].append(zpe * per_mol / 1000)
            continue
        x = np.clip(e / (KB_EV * T), 1e-12, 700)
        n = 1 / np.expm1(x)                              # Bose occupation
        U = np.sum(g * dw * e * (0.5 + n))
        F = np.sum(g * dw * (e / 2 + KB_EV * T * np.log1p(-np.exp(-x))))
        Cv = np.sum(g * dw * KB_EV * x ** 2 * np.exp(-x) / (-np.expm1(-x)) ** 2)  # eV/K per cell
        S = (U - F) / T
        out["Cv"].append(Cv * per_mol)                  # J/(mol K)
        out["S"].append(S * per_mol)                    # J/(mol K)
        out["F"].append(F * per_mol / 1000)             # kJ/mol
        out["U"].append(U * per_mol / 1000)             # kJ/mol
    out["zpe_ev_atom"] = zpe / natoms
    out["imaginary_fraction"] = float(neg)
    out["method"] = "from the phonon DOS (broadened; entropy slightly overestimated)"
    out["units"] = {"Cv": "J/(mol·K)", "S": "J/(mol·K)", "F": "kJ/mol", "U": "kJ/mol", "per": "mole of atoms"}
    return out


def thermodynamics_modes(freq_cm: np.ndarray, natoms: int, temps: np.ndarray | None = None) -> dict:
    """Harmonic thermodynamics per mole of atoms by direct summation over all modes of a uniform q-mesh.

    More accurate than integrating a broadened DOS, which shifts weight to low frequencies and overestimates
    the entropy and the low-temperature heat capacity. freq_cm: (nq, 3*natoms).
    """
    if temps is None:
        temps = np.arange(0.0, 1001.0, 10.0)
    f = np.asarray(freq_cm, float)
    nq = f.shape[0]
    allf = f.ravel()
    imag = float(np.mean(allf < IMAG_THRESHOLD_CM1))
    e = allf[allf > 1e-3] * CM1_EV                       # eV; acoustic zeros at Gamma and imaginary modes dropped
    per_mol = EV_JMOL / natoms / nq                       # eV summed over the mesh -> J per mole of atoms
    zpe = np.sum(e / 2)
    out = {"T": temps.tolist(), "Cv": [], "S": [], "F": [], "U": []}
    for T in temps:
        if T < 1e-8:
            out["Cv"].append(0.0)
            out["S"].append(0.0)
            out["F"].append(zpe * per_mol / 1000)
            out["U"].append(zpe * per_mol / 1000)
            continue
        x = np.clip(e / (KB_EV * T), 1e-12, 700)
        n = 1 / np.expm1(x)
        U = np.sum(e * (0.5 + n))
        F = np.sum(e / 2 + KB_EV * T * np.log1p(-np.exp(-x)))
        Cv = np.sum(KB_EV * x ** 2 * np.exp(-x) / (-np.expm1(-x)) ** 2)
        out["Cv"].append(Cv * per_mol)
        out["S"].append((U - F) / T * per_mol)
        out["F"].append(F * per_mol / 1000)
        out["U"].append(U * per_mol / 1000)
    out["zpe_ev_atom"] = zpe / natoms / nq
    out["imaginary_fraction"] = imag
    out["method"] = f"direct sum over a {round(nq ** (1 / 3))}³ q-mesh"
    out["units"] = {"Cv": "J/(mol·K)", "S": "J/(mol·K)", "F": "kJ/mol", "U": "kJ/mol", "per": "mole of atoms"}
    return out


def value_at(thermo: dict, key: str, T: float) -> float | None:
    Ts = np.array(thermo["T"])
    if not len(Ts):
        return None
    return float(np.interp(T, Ts, np.array(thermo[key])))
