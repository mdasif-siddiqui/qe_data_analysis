"""Optical spectra (UV-Vis) in the independent-particle approximation, straight from QE's HDF5 wavefunctions.

epsilon.x of this QE build gives k-grid dependent, several-times-too-large spectra, so the dielectric function
is computed here instead, with the same physics (RPA without local fields, no [V_NL, r] commutator):

    eps2(w) = 4 pi^2 / (Omega w^2) * sum_k w_k sum_{n<m} (f_n - f_m) |<m|p|n>|^2 / 3 * delta(E_m - E_n - w)

averaged over the three directions. The trace is invariant under the crystal's symmetry operations, so the
symmetry-reduced k-points of a normal NSCF run give the exact orientational average. eps1 follows from
Kramers-Kronig. Everything in Hartree atomic units internally.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

HA_EV = 27.211386245988
_trapz = getattr(np, "trapezoid", None) or np.trapz
HC_EV_NM = 1239.84198


def read_xml(save: Path) -> dict:
    root = ET.parse(save / "data-file-schema.xml").getroot()
    out = root.find("output")
    bs = out.find("band_structure")
    lsda = bs.findtext("lsda", "false").strip() == "true"
    noncolin = bs.findtext("noncolin", "false").strip() == "true"
    nelec = float(bs.findtext("nelec"))
    cell = out.find("atomic_structure/cell")
    a = np.array([[float(x) for x in cell.findtext(t).split()] for t in ("a1", "a2", "a3")])
    vol = abs(np.linalg.det(a))
    ks = bs.findall("ks_energies")
    w = np.array([float(k.find("k_point").attrib["weight"]) for k in ks])
    alat = float(out.find("atomic_structure").attrib.get("alat", 1.0))
    kc = np.array([[float(x) for x in k.findtext("k_point").split()] for k in ks]) * 2 * np.pi / alat   # 1/bohr
    eig = np.array([[float(x) for x in k.findtext("eigenvalues").split()] for k in ks])        # Ha
    occ = np.array([[float(x) for x in k.findtext("occupations").split()] for k in ks])        # 0..1
    kind = (bs.findtext("occupations_kind") or "").strip().lower()
    ef = bs.findtext("fermi_energy")
    return {"lsda": lsda, "noncolin": noncolin, "nelec": nelec, "volume": vol, "w": w / w.sum(),
            "eig": eig, "occ": occ, "k": kc, "smearing": kind == "smearing",
            "fermi": float(ef) if ef is not None and kind == "smearing" else None}


def read_wfc(path: Path) -> tuple[np.ndarray, np.ndarray]:
    import h5py
    with h5py.File(path, "r") as f:
        mill = f["MillerIndices"][:]
        bg = np.array([f["MillerIndices"].attrs[f"bg{i}"] for i in (1, 2, 3)])
        ev = f["evc"][:]
    c = ev[:, 0::2] + 1j * ev[:, 1::2]                    # (nbnd, npw*npol)
    G = mill @ bg                                          # Cartesian, 1/bohr
    return c, G


def _channels(save: Path, x: dict):
    """(wfc file prefix, slice of bands in the XML arrays, spin degeneracy) per spin channel."""
    nb = x["eig"].shape[1]
    if x["lsda"]:
        half = nb // 2
        return [("wfcup", slice(0, half), 1.0), ("wfcdw", slice(half, nb), 1.0)]
    return [("wfc", slice(0, nb), 1.0 if x["noncolin"] else 2.0)]


def dielectric(save: str | Path, emax_ev: float = 20.0, sigma_ev: float = 0.1, de_ev: float = 0.01,
               drude_width_ev: float = 0.3) -> dict:
    """Interband part as described in the module docstring. For metals (smearing, Fermi level in a band) the
    intraband (Drude) plasma frequency is added:

        wp^2 = 4 pi / Omega * sum_k w_k g sum_n |<n|k+G|n>|^2 / 3 * delta(E_n - E_F)

    (free electrons: exactly 4 pi N / Omega). The delta function is a Gaussian of width drude_width_ev,
    which must be wide compared with the k-point spacing near E_F."""
    save = Path(save)
    x = read_xml(save)
    E = np.arange(0.0, emax_ev + de_ev / 2, de_ev) / HA_EV           # Ha
    sig = sigma_ev / HA_EV
    eps2 = np.zeros_like(E)
    wd = drude_width_ev / HA_EV
    wp2_intra = 0.0
    for prefix, sl, g in _channels(save, x):
        for ik in range(x["eig"].shape[0]):
            c, G = read_wfc(save / f"{prefix}{ik + 1}.hdf5")
            npw = G.shape[0]
            if c.shape[1] != npw:                                     # spinors: two components share G
                c = c.reshape(c.shape[0], -1, npw).transpose(0, 2, 1).reshape(c.shape[0], npw, -1)
                p = [np.einsum("mgs,ngs->mn", c.conj(), c * G[None, :, a, None]) for a in range(3)]
                dens = (np.abs(c) ** 2).sum(axis=2)                       # (nbnd, npw)
            else:
                p = [(c.conj() * G[:, a]) @ c.T for a in range(3)]   # <m|p_a|n> (off-diagonal; k+G -> G)
                dens = np.abs(c) ** 2
            p2 = (np.abs(p[0]) ** 2 + np.abs(p[1]) ** 2 + np.abs(p[2]) ** 2) / 3
            e = x["eig"][ik, sl]
            f = x["occ"][ik, sl]
            if x["fermi"] is not None:
                # intraband: diagonal velocity <n|k+G|n> (k matters here), states within ~wd of E_F
                pnn = dens @ (G + x["k"][ik])                             # (nbnd, 3)
                near = np.exp(-0.5 * ((e - x["fermi"]) / wd) ** 2) / (wd * np.sqrt(2 * np.pi))
                wp2_intra += x["w"][ik] * g * float(np.sum(near * (pnn ** 2).sum(axis=1) / 3))
            dE = e[None, :] - e[:, None]                               # E_m - E_n  (rows m? -> use [n, m])
            df = f[:, None] - f[None, :]                               # f_n - f_m
            mask = (dE > 1e-6) & (df > 1e-8)
            n_idx, m_idx = np.nonzero(mask)
            if not len(n_idx):
                continue
            wt = x["w"][ik] * g * df[n_idx, m_idx] * p2[m_idx, n_idx]
            om = dE[n_idx, m_idx]
            # Gaussian-broadened delta, eps2 ~ 1/w^2 * sum
            gauss = np.exp(-0.5 * ((E[:, None] - om[None, :]) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
            eps2 += gauss @ wt
    with np.errstate(divide="ignore", invalid="ignore"):
        eps2 = np.where(E > 1e-9, 4 * np.pi ** 2 / x["volume"] * eps2 / E ** 2, 0.0)
    eps1 = kramers_kronig(E, eps2)
    # f-sum rule: int w eps2 dw = (pi/2) wp^2, wp^2 = 4 pi N/Omega. Ratio < 1 means too few empty bands.
    wp2 = 4 * np.pi * x["nelec"] / x["volume"]
    wp2_intra *= 4 * np.pi / x["volume"]
    # the intraband (Drude) part carries exactly (pi/2) wp_intra^2 of the sum rule
    fsum = (_trapz(E * eps2, E) + np.pi / 2 * wp2_intra) / (np.pi / 2 * wp2)
    return {"E_ev": E * HA_EV, "eps1": eps1, "eps2": eps2, "fsum_ratio": float(fsum),
            "plasma_ev": float(np.sqrt(wp2) * HA_EV), "nelec": x["nelec"], "nbnd": int(x["eig"].shape[1]),
            "nk": int(x["eig"].shape[0]), "metal": bool(x["fermi"] is not None and wp2_intra > (0.5 / HA_EV) ** 2),
            "drude_plasma_ev": float(np.sqrt(wp2_intra) * HA_EV), "drude_width_ev": drude_width_ev}


def drude(E_ev: np.ndarray, wp_ev: float, gamma_ev: float) -> tuple[np.ndarray, np.ndarray]:
    """Free-carrier term eps_D = -wp^2 / (w^2 + i gamma w)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        d = E_ev ** 2 + gamma_ev ** 2
        e1 = -wp_ev ** 2 / d
        e2 = np.where(E_ev > 0, wp_ev ** 2 * gamma_ev / (E_ev * d), 0.0)
    return e1, e2


def kramers_kronig(E: np.ndarray, eps2: np.ndarray) -> np.ndarray:
    """eps1(w) = 1 + 2/pi P int w' eps2(w') / (w'^2 - w^2) dw' (singular point left out, uniform grid)."""
    dE = E[1] - E[0]
    eps1 = np.ones_like(E)
    num = E * eps2
    for i, w in enumerate(E):
        d = E ** 2 - w ** 2
        d[i] = np.inf
        eps1[i] += 2 / np.pi * np.sum(num / d) * dE
    return eps1


def optical_constants(E_ev: np.ndarray, eps1: np.ndarray, eps2: np.ndarray) -> dict:
    mod = np.sqrt(eps1 ** 2 + eps2 ** 2)
    n = np.sqrt(np.clip((mod + eps1) / 2, 0, None))
    k = np.sqrt(np.clip((mod - eps1) / 2, 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        lam_cm = np.where(E_ev > 0, HC_EV_NM / E_ev * 1e-7, np.inf)
        alpha = np.where(E_ev > 0, 4 * np.pi * k / lam_cm, 0.0)         # 1/cm
    R = ((n - 1) ** 2 + k ** 2) / ((n + 1) ** 2 + k ** 2)
    return {"n": n, "k": k, "alpha_cm": alpha, "R": R}


def scissor(E: np.ndarray, eps2: np.ndarray, shift_ev: float) -> np.ndarray:
    """Rigid shift of the absorption edge: eps2_s(w) = (w - D)/w * eps2(w - D) (keeps the f-sum rule)."""
    if not shift_ev:
        return eps2
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.interp(E - shift_ev, E, eps2, left=0.0) * np.where(E > 0, (E - shift_ev) / E, 0.0)
    return np.clip(out, 0, None)


def compute(save: str | Path, outdir: str | Path, emax_ev: float = 20.0, sigma_ev: float = 0.1,
            shift_ev: float = 0.0) -> dict:
    d = dielectric(save, emax_ev=emax_ev, sigma_ev=sigma_ev)
    E = d["E_ev"]
    eps2 = scissor(E, d["eps2"], shift_ev)
    eps1 = kramers_kronig(E / HA_EV, eps2) if shift_ev else d["eps1"]
    eps_inter_static = float(eps1[0])
    if d["metal"]:
        # free carriers (intraband); damping = the broadening, so the Drude peak has the same resolution
        d1, d2 = drude(E, d["drude_plasma_ev"], sigma_ev)
        eps1, eps2 = eps1 + d1, eps2 + d2
        eps1[0], eps2[0] = eps1[1], eps2[1]                               # w = 0 is singular
    oc = optical_constants(E, eps1, eps2)
    outdir = Path(outdir)
    np.savetxt(outdir / "optics.dat", np.column_stack([E, eps1, eps2, oc["n"], oc["k"], oc["alpha_cm"], oc["R"]]),
               header="E(eV) eps1 eps2 n k alpha(1/cm) R   [independent-particle, orientation average]",
               fmt="%.6e")
    meta = {k: d[k] for k in ("fsum_ratio", "plasma_ev", "nelec", "nbnd", "nk", "metal", "drude_plasma_ev", "drude_width_ev")}
    meta.update(sigma_ev=sigma_ev, shift_ev=shift_ev, emax_ev=emax_ev,
                eps_static=None if d["metal"] else eps_inter_static, eps_interband_static=eps_inter_static)
    if d["metal"]:
        meta["drude_gamma_ev"] = sigma_ev
    return meta
