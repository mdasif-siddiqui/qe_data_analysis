"""Optical spectra (UV-Vis) in the independent-particle approximation, straight from QE's HDF5 wavefunctions.

epsilon.x of this QE build gives k-grid dependent, several-times-too-large spectra, so the dielectric function
is computed here instead (RPA without local fields):

    eps2(w) = 4 pi^2 / (Omega w^2) * sum_k w_k sum_{n<m} (f_n - f_m) |<m|v|n>|^2 / 3 * delta(E_m - E_n - w)

with the velocity v = dH_k/dk = p + i[V_NL, r]. The nonlocal part dV_NL(k)/dk is built from the
pseudopotentials' beta projectors (norm-conserving only) by central differences in k; leaving it out
overestimates eps by 10-20%. The trace is invariant under the crystal's symmetry operations, so the
symmetry-reduced k-points of a normal NSCF run give the exact orientational average. eps1 follows from
Kramers-Kronig. Everything in Hartree atomic units internally.
"""
from __future__ import annotations

import re
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


# ----------------------------------------------------------------------------- nonlocal pseudopotential
def _floats(s: str) -> np.ndarray:
    return np.array([float(t) for t in s.replace(",", " ").split()])


def _simpson(f: np.ndarray, rab: np.ndarray) -> float:
    """QE-style Simpson integration on the radial mesh (rab = dr/di)."""
    n = len(f)
    w = np.ones(n)
    w[1:-1:2], w[2:-1:2] = 4.0, 2.0
    if n % 2 == 0:                       # even number of points: last interval by the trapezoid rule
        w[-1] = 1.0
    return float(np.sum(w * f * rab) / 3.0)


def read_upf_nonlocal(path: str | Path) -> dict | None:
    """beta projectors (r*beta(r)), their l, D_ij (Ry) and the radial mesh of a norm-conserving UPF (v1 or v2).

    None for ultrasoft / PAW files (their [S, r] terms are not implemented)."""
    txt = Path(path).read_text(errors="replace")
    head = re.search(r"<PP_HEADER(.*?)>(.*?)(?:</PP_HEADER>|(?=<PP_MESH))", txt, re.S)
    htxt = (head.group(1) + head.group(2)) if head else ""
    if re.search(r'pseudo_type\s*=\s*"\s*(US|PAW|USPP)', htxt, re.I) or re.search(r"\b(US|PAW)\s+Ultrasoft|^\s*(US|PAW)\s", htxt, re.M):
        return None
    if re.search(r'is_ultrasoft\s*=\s*"\s*(t|true|\.true\.)', htxt, re.I) or re.search(r'is_paw\s*=\s*"\s*(t|true|\.true\.)', htxt, re.I):
        return None
    r = _floats(re.search(r"<PP_R[^>]*>(.*?)</PP_R>", txt, re.S).group(1))
    rab = _floats(re.search(r"<PP_RAB[^>]*>(.*?)</PP_RAB>", txt, re.S).group(1))
    betas = []
    v2 = list(re.finditer(r"<PP_BETA\.(\d+)([^>]*)>(.*?)</PP_BETA\.\1>", txt, re.S))
    if v2:
        for m in v2:
            l = int(re.search(r'angular_momentum\s*=\s*"\s*(\d+)', m.group(2)).group(1))
            kk = re.search(r'cutoff_radius_index\s*=\s*"\s*(\d+)', m.group(2))
            vals = _floats(m.group(3))
            kk = int(kk.group(1)) if kk and int(kk.group(1)) > 0 else len(vals)
            betas.append((l, vals[:kk]))
        dij = _floats(re.search(r"<PP_DIJ[^>]*>(.*?)</PP_DIJ>", txt, re.S).group(1))
        nb = len(betas)
        D = dij[: nb * nb].reshape(nb, nb)
    else:                                                   # UPF v1
        for m in re.finditer(r"<PP_BETA>(.*?)</PP_BETA>", txt, re.S):
            lines = m.group(1).strip().splitlines()
            l = int(lines[0].split()[1])
            kk = int(lines[1].split()[0])
            vals = _floats(" ".join(lines[2:]))[:kk]
            betas.append((l, vals))
        nb = len(betas)
        D = np.zeros((nb, nb))
        blk = re.search(r"<PP_DIJ>(.*?)</PP_DIJ>", txt, re.S)
        if blk:
            lines = blk.group(1).strip().splitlines()
            for ln in lines[1:]:
                p = ln.split()
                if len(p) >= 3:
                    i, j = int(p[0]) - 1, int(p[1]) - 1
                    D[i, j] = D[j, i] = float(p[2].lower().replace("d", "e"))
    return {"r": r, "rab": rab, "betas": betas, "D": D}


def _ylm(l: int, u: np.ndarray) -> np.ndarray:
    """Real orthonormal spherical harmonics, (2l+1, N), of unit vectors u (N, 3)."""
    x, y, z = u[:, 0], u[:, 1], u[:, 2]
    pi = np.pi
    if l == 0:
        return np.full((1, len(u)), 1 / np.sqrt(4 * pi))
    if l == 1:
        return np.sqrt(3 / (4 * pi)) * np.array([x, y, z])
    if l == 2:
        c = np.sqrt(15 / (4 * pi))
        return np.array([c * x * y, c * y * z, np.sqrt(5 / (16 * pi)) * (3 * z * z - 1), c * x * z, c / 2 * (x * x - y * y)])
    if l == 3:
        return np.array([np.sqrt(35 / (32 * pi)) * (3 * x * x - y * y) * y, np.sqrt(105 / (4 * pi)) * x * y * z,
                         np.sqrt(21 / (32 * pi)) * y * (5 * z * z - 1), np.sqrt(7 / (16 * pi)) * z * (5 * z * z - 3),
                         np.sqrt(21 / (32 * pi)) * x * (5 * z * z - 1), np.sqrt(105 / (16 * pi)) * z * (x * x - y * y),
                         np.sqrt(35 / (32 * pi)) * (x * x - 3 * y * y) * x])
    raise ValueError(f"l = {l} projectors are not supported")


class NonlocalPotential:
    """V_NL(k) = sum_I sum_ij |beta_i^I> D_ij <beta_j^I| in the plane-wave basis of one k-point (Hartree)."""

    def __init__(self, species: dict[str, dict], atoms: list[tuple[str, np.ndarray]], volume: float, qmax: float):
        from scipy.interpolate import CubicSpline
        from scipy.special import spherical_jn
        self.vol = volume
        self.atoms = atoms
        self.tab = {}
        q = np.arange(0.0, qmax + 0.5, 0.01)
        for name, upf in species.items():
            r, rab = upf["r"], upf["rab"]
            splines = []
            for l, rb in upf["betas"]:
                n = len(rb)
                f = np.array([_simpson(r[:n] * rb * spherical_jn(l, qq * r[:n]), rab[:n]) for qq in q])
                splines.append((l, CubicSpline(q, f)))
            self.tab[name] = (splines, upf["D"] / 2.0)       # D_ij: Ry -> Ha

    def projectors(self, kG: np.ndarray, G: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(B, D): B (npw, nproj) with B_p(k+G) = 4 pi / sqrt(Omega) f(|k+G|) Y_lm e^{-i G.tau}; D block-diagonal."""
        qn = np.linalg.norm(kG, axis=1)
        u = kG / np.where(qn > 1e-12, qn, 1.0)[:, None]
        u[qn <= 1e-12] = [0.0, 0.0, 1.0]
        cols, blocks = [], []
        for name, tau in self.atoms:
            splines, D = self.tab[name]
            phase = np.exp(-1j * (G @ tau))
            idx = []
            for i, (l, sp) in enumerate(splines):
                y = _ylm(l, u)
                fq = sp(qn)
                for m in range(2 * l + 1):
                    cols.append(4 * np.pi / np.sqrt(self.vol) * fq * y[m] * phase)
                    idx.append((i, l, m))
            nb = len(idx)
            Dblk = np.zeros((nb, nb))
            for a, (i, li, mi) in enumerate(idx):
                for b, (j, lj, mj) in enumerate(idx):
                    if li == lj and mi == mj:
                        Dblk[a, b] = D[i, j]
            blocks.append(Dblk)
        B = np.array(cols).T if cols else np.zeros((len(G), 0))
        n = sum(b.shape[0] for b in blocks)
        Dall = np.zeros((n, n))
        o = 0
        for b in blocks:
            Dall[o:o + b.shape[0], o:o + b.shape[0]] = b
            o += b.shape[0]
        return B, Dall

    def velocity(self, c: np.ndarray, k: np.ndarray, G: np.ndarray, dk: float = 1e-4) -> list[np.ndarray]:
        """<m| dV_NL/dk_a |n> for a = x, y, z (central differences), c: (nbnd, npw)."""
        out = []
        for a in range(3):
            e = np.zeros(3)
            e[a] = dk
            V = []
            for s in (+1, -1):
                B, D = self.projectors(k + s * e + G, G)
                A = B.conj().T @ c.T                           # (nproj, nbnd): <beta_p|n>
                V.append(A.conj().T @ D @ A)
            out.append((V[0] - V[1]) / (2 * dk))
        return out


def nonlocal_from_save(save: Path, extra_dirs: list[Path] | None = None) -> tuple[NonlocalPotential | None, str]:
    """Build the nonlocal potential from the save directory (positions, species) and the UPF files."""
    root = ET.parse(save / "data-file-schema.xml").getroot()
    out = root.find("output")
    if (out.findtext("band_structure/noncolin") or "false").strip() == "true":
        return None, "noncollinear run: the nonlocal commutator is not included"
    a = np.array([[float(x) for x in out.findtext(f"atomic_structure/cell/{t}").split()] for t in ("a1", "a2", "a3")])
    vol = abs(np.linalg.det(a))
    species = {}
    for sp in out.findall("atomic_species/species"):
        name, pf = sp.attrib["name"], sp.findtext("pseudo_file").strip()
        path = next((dd / pf for dd in [save] + list(extra_dirs or []) if (dd / pf).exists()), None)
        if path is None:
            return None, f"pseudopotential {pf} not found: the nonlocal commutator is not included"
        upf = read_upf_nonlocal(path)
        if upf is None:
            return None, f"{pf} is ultrasoft/PAW: the nonlocal commutator is not included"
        species[name] = upf
    atoms = [(at.attrib["name"], np.array([float(x) for x in at.text.split()]))
             for at in out.findall("atomic_structure/atomic_positions/atom")]
    ecut = float(out.findtext("basis_set/ecutwfc"))                    # Ha
    return NonlocalPotential(species, atoms, vol, qmax=np.sqrt(2 * ecut) + 1.0), "included"


def _channels(save: Path, x: dict):
    """(wfc file prefix, slice of bands in the XML arrays, spin degeneracy) per spin channel."""
    nb = x["eig"].shape[1]
    if x["lsda"]:
        half = nb // 2
        return [("wfcup", slice(0, half), 1.0), ("wfcdw", slice(half, nb), 1.0)]
    return [("wfc", slice(0, nb), 1.0 if x["noncolin"] else 2.0)]


def dielectric(save: str | Path, emax_ev: float = 20.0, sigma_ev: float = 0.1, de_ev: float = 0.01,
               drude_width_ev: float = 0.3, nonlocal_: bool = True, pseudo_dirs: list | None = None) -> dict:
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
    fsum_p = 0.0
    vnl, nl_status = (None, "switched off")
    if nonlocal_:
        try:
            vnl, nl_status = nonlocal_from_save(save, [Path(p) for p in (pseudo_dirs or [])])
        except Exception as e:                                        # never lose the spectrum over this term
            vnl, nl_status = None, f"could not be built ({e}): not included"
    for prefix, sl, g in _channels(save, x):
        for ik in range(x["eig"].shape[0]):
            c, G = read_wfc(save / f"{prefix}{ik + 1}.hdf5")
            npw = G.shape[0]
            dV = None
            if c.shape[1] != npw:                                     # spinors: two components share G
                c = c.reshape(c.shape[0], -1, npw).transpose(0, 2, 1).reshape(c.shape[0], npw, -1)
                p = [np.einsum("mgs,ngs->mn", c.conj(), c * G[None, :, a, None]) for a in range(3)]
                dens = (np.abs(c) ** 2).sum(axis=2)                       # (nbnd, npw)
            else:
                p = [(c.conj() * G[:, a]) @ c.T for a in range(3)]   # <m|p_a|n> (off-diagonal; k+G -> G)
                dens = np.abs(c) ** 2
            p2_bare = (np.abs(p[0]) ** 2 + np.abs(p[1]) ** 2 + np.abs(p[2]) ** 2) / 3
            if dV is None and vnl is not None and c.ndim == 2:        # v = p + dV_NL/dk
                dV = vnl.velocity(c, x["k"][ik], G)
                p = [p[a] + dV[a] for a in range(3)]
            p2 = (np.abs(p[0]) ** 2 + np.abs(p[1]) ** 2 + np.abs(p[2]) ** 2) / 3
            e = x["eig"][ik, sl]
            f = x["occ"][ik, sl]
            if x["fermi"] is not None:
                # intraband: diagonal velocity <n|k+G|n> (+ <n|dV_NL/dk|n>), states within ~wd of E_F
                pnn = dens @ (G + x["k"][ik])                             # (nbnd, 3)
                if dV is not None:
                    pnn = pnn + np.array([np.real(np.diag(dV[a])) for a in range(3)]).T
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
            # f-sum (band-count check) with the bare momentum: the TRK sum rule is exact only for p
            fsum_p += float(np.sum(x["w"][ik] * g * df[n_idx, m_idx] * p2_bare[m_idx, n_idx] / om))
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
    # int w eps2 dw = 4 pi^2/Omega sum wt/om (exact for the Gaussian-broadened deltas)
    fsum = (4 * np.pi ** 2 / x["volume"] * fsum_p + np.pi / 2 * wp2_intra) / (np.pi / 2 * wp2)
    return {"E_ev": E * HA_EV, "eps1": eps1, "eps2": eps2, "fsum_ratio": float(fsum),
            "plasma_ev": float(np.sqrt(wp2) * HA_EV), "nelec": x["nelec"], "nbnd": int(x["eig"].shape[1]),
            "nk": int(x["eig"].shape[0]), "metal": bool(x["fermi"] is not None and wp2_intra > (0.5 / HA_EV) ** 2),
            "drude_plasma_ev": float(np.sqrt(wp2_intra) * HA_EV), "drude_width_ev": drude_width_ev,
            "nonlocal": nl_status}


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
            shift_ev: float = 0.0, pseudo_dirs: list | None = None) -> dict:
    d = dielectric(save, emax_ev=emax_ev, sigma_ev=sigma_ev,
                   pseudo_dirs=list(pseudo_dirs or []) + [Path(outdir) / "pseudos", Path(outdir).parent / "pseudos"])
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
    meta = {k: d[k] for k in ("fsum_ratio", "plasma_ev", "nelec", "nbnd", "nk", "metal", "drude_plasma_ev", "drude_width_ev",
                              "nonlocal")}
    meta.update(sigma_ev=sigma_ev, shift_ev=shift_ev, emax_ev=emax_ev,
                eps_static=None if d["metal"] else eps_inter_static, eps_interband_static=eps_inter_static)
    if d["metal"]:
        meta["drude_gamma_ev"] = sigma_ev
    return meta
