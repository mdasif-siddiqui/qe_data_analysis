"""Convert a CIF file into a ready-to-run pw.x scf input."""
from __future__ import annotations

import io
import math
import re
import urllib.error
import urllib.request

import numpy as np

from . import pseudos

# Pseudopotential file names tried in order on the QE library, per functional.
# "lda-nc" (norm-conserving LDA) is what Raman (ph.x lraman) and UV-Vis (epsilon.x) need.
CANDIDATES = {
    "pbe": ["{el}.pbe-n-rrkjus_psl.1.0.0.UPF", "{el}.pbe-dn-rrkjus_psl.1.0.0.UPF", "{el}.pbe-spn-rrkjus_psl.1.0.0.UPF",
            "{el}.pbe-spdn-rrkjus_psl.1.0.0.UPF", "{el}.pbe-n-kjpaw_psl.1.0.0.UPF", "{el}.pbe-dn-kjpaw_psl.1.0.0.UPF",
            "{el}.pbe-spn-kjpaw_psl.1.0.0.UPF", "{el}.pbe-spdn-kjpaw_psl.1.0.0.UPF",
            "{el}.pbe-rrkjus_psl.1.0.0.UPF", "{el}.pbe-kjpaw_psl.1.0.0.UPF", "{el}.pbe-rrkjus.UPF"],  # H, He: no valence tag
    "pbe-paw": ["{el}.pbe-n-kjpaw_psl.1.0.0.UPF", "{el}.pbe-dn-kjpaw_psl.1.0.0.UPF", "{el}.pbe-spn-kjpaw_psl.1.0.0.UPF",
                "{el}.pbe-spdn-kjpaw_psl.1.0.0.UPF", "{el}.pbe-kjpaw_psl.1.0.0.UPF", "{el}.pbe-nl-kjpaw_psl.1.0.0.UPF"],
    "lda-nc": ["{el}.pz-vbc.UPF", "{el}.pz-bhs.UPF", "{el}.pz-mt.UPF", "{el}.pz-n-vbc.UPF", "{el}.pz-hgh.UPF"],
}
FUNCTIONALS = {"pbe": "PBE (ultrasoft/PAW, PSlibrary)", "pbe-paw": "PBE PAW (PSlibrary; GIPAW-ready for NMR)", "lda-nc": "LDA, norm-conserving (needed for Raman and UV-Vis)"}
_exists_cache: dict[str, bool] = {}


def _exists_online(name: str) -> bool:
    if name in _exists_cache:
        return _exists_cache[name]
    req = urllib.request.Request(pseudos.QE_URL.format(name=name), method="HEAD", headers={"User-Agent": "qe-studio"})
    try:
        with urllib.request.urlopen(req, timeout=15) as r:
            ok = r.status == 200
    except urllib.error.HTTPError:
        ok = False
    except (urllib.error.URLError, TimeoutError, OSError):
        return False  # offline: do not cache
    _exists_cache[name] = ok
    return ok


def choose_pseudo(el: str, functional: str) -> str | None:
    """First candidate that is already in the library, else the first that exists online."""
    names = [c.format(el=el) for c in CANDIDATES[functional]]
    for n in names:
        if pseudos.find(n):
            return n
    for n in names:
        if _exists_online(n):
            return n
    return None


def suggested_cutoffs(path) -> tuple[float | None, float | None]:
    """PSlibrary files state 'Suggested minimum cutoff for wavefunctions/charge density' in PP_INFO."""
    txt = path.read_text(errors="replace")[:20000]
    w = re.search(r"Suggested minimum cutoff for wavefunctions:\s*([\d.]+)", txt)
    r = re.search(r"Suggested minimum cutoff for charge density:\s*([\d.]+)", txt)
    return (float(w.group(1)) if w else None, float(r.group(1)) if r else None)


def pseudo_kind(path) -> str:
    head = path.read_text(errors="replace")[:20000]
    m = re.search(r'pseudo_type\s*=\s*"(\w+)"', head)
    if m:
        return m.group(1).upper()
    if re.search(r"\bUS\b|ultrasoft", head, re.I):
        return "US"
    return "NC"


def cif_to_scf(cif_text: str, functional: str = "pbe", primitive: bool = True, occupations: str = "smearing",
               kspacing: float = 0.25, prefix: str | None = None, download: bool = True) -> dict:
    """Returns {"text": scf input, "info": {...}, "warnings": [...]}."""
    from ase.io import read
    from ase import Atoms
    warnings: list[str] = []
    if functional not in CANDIDATES:
        raise ValueError(f"Unknown functional {functional!r}")
    if re.search(r"_atom_site_occupancy", cif_text):
        occ = re.findall(r"\b0\.\d+\b", cif_text.split("_atom_site_occupancy", 1)[1][:4000])
        if occ:
            warnings.append("The CIF lists partial site occupancies (disorder). Each site is treated as fully occupied; "
                            "an ordered supercell is needed for a realistic calculation.")
    try:
        atoms = read(io.StringIO(cif_text), format="cif")
    except Exception as e:
        raise ValueError(f"Could not read the CIF file: {e}")
    if isinstance(atoms, list):
        atoms = atoms[0]
    info: dict = {"formula_cif": atoms.get_chemical_formula(mode="hill"), "natoms_cif": len(atoms)}
    try:
        import spglib
        cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
        ds = spglib.get_symmetry_dataset(cell, symprec=1e-3)
        if ds is not None:
            info["spacegroup"] = f"{ds.international} (No. {ds.number})"
        if primitive:
            lat, pos, nums = spglib.standardize_cell(cell, to_primitive=True, no_idealize=False, symprec=1e-3)
            if len(nums) < len(atoms):
                atoms = Atoms(numbers=nums, cell=lat, scaled_positions=pos, pbc=True)
    except ImportError:
        warnings.append("spglib is not installed: the cell was not reduced and the space group is unknown.")
    info["natoms"] = len(atoms)
    info["formula"] = atoms.get_chemical_formula(mode="hill")
    elements = list(dict.fromkeys(atoms.get_chemical_symbols()))
    from ase.data import atomic_masses, atomic_numbers

    # pseudopotentials and cutoffs
    ps, ecw, ecr, kinds = {}, [], [], set()
    for el in elements:
        name = choose_pseudo(el, functional)
        ps[el] = name
        if not name:
            warnings.append(f"No {FUNCTIONALS[functional]} pseudopotential for {el} found in the QE library: "
                            f"upload a .UPF file for {el} and edit ATOMIC_SPECIES.")
            continue
        path = pseudos.find(name)
        if path is None and download:
            try:
                path = pseudos.download(name)
            except Exception as e:
                warnings.append(f"Could not download {name}: {e}")
        if path:
            w, r = suggested_cutoffs(path)
            kinds.add(pseudo_kind(path))
            if w:
                ecw.append(w)
            if r:
                ecr.append(r)
            if "hgh" in name:
                ecw.append(80.0)
    nc_only = kinds <= {"NC", "SL"}
    ecutwfc = max(ecw) if ecw else (50.0 if nc_only else 45.0)
    ecutwfc = float(math.ceil(ecutwfc / 5) * 5)
    ecutrho = max(ecr) if ecr else (4 * ecutwfc if nc_only else 8 * ecutwfc)
    ecutrho = float(max(math.ceil(ecutrho / 10) * 10, (4 if nc_only else 8) * ecutwfc))
    if not ecw:
        warnings.append(f"The pseudopotentials do not state a recommended cutoff; ecutwfc = {ecutwfc:g} Ry is a guess. "
                        "Run the cutoff convergence test.")
    info.update(ecutwfc=ecutwfc, ecutrho=ecutrho, pseudos=ps, pseudo_types=sorted(kinds), functional=FUNCTIONALS[functional])

    # k-points from a target spacing (2*pi/a convention, 1/Angstrom)
    rec = np.linalg.norm(atoms.cell.reciprocal(), axis=1) * 2 * math.pi
    grid = [max(1, int(math.ceil(b / kspacing))) for b in rec]
    info["kgrid"] = grid

    prefix = prefix or re.sub(r"[^\w]", "", info["formula"]) or "pwscf"
    lines = ["&CONTROL", "  calculation = 'scf'", f"  prefix = '{prefix}'", "  outdir = './tmp/'",
             "  pseudo_dir = './pseudos/'", "  tprnfor = .true.", "  tstress = .true.", "/", "&SYSTEM", "  ibrav = 0",
             f"  nat = {len(atoms)}", f"  ntyp = {len(elements)}", f"  ecutwfc = {ecutwfc:g}", f"  ecutrho = {ecutrho:g}"]
    if occupations == "smearing":
        lines += ["  occupations = 'smearing'", "  smearing = 'mv'", "  degauss = 0.01"]
    lines += ["/", "&ELECTRONS", "  conv_thr = 1.0d-9", "  mixing_beta = 0.5", "/", "ATOMIC_SPECIES"]
    for el in elements:
        mass = atomic_masses[atomic_numbers[el]]
        lines.append(f"  {el:2s} {mass:9.4f} {ps[el] or el + '.UPF'}")
    lines.append("CELL_PARAMETERS angstrom")
    for v in atoms.cell.array:
        lines.append("  " + " ".join(f"{x:14.9f}" for x in v))
    lines.append("ATOMIC_POSITIONS crystal")
    for s, p in zip(atoms.get_chemical_symbols(), atoms.get_scaled_positions(wrap=True)):
        lines.append(f"  {s:2s} " + " ".join(f"{x:13.9f}" for x in p))
    lines += ["K_POINTS automatic", "  {} {} {} 0 0 0".format(*grid)]
    if occupations == "fixed":
        warnings.append("Fixed occupations: correct only for insulators and semiconductors. Use smearing for metals.")
    if len(atoms) > 40:
        warnings.append(f"{len(atoms)} atoms: this will be slow on a laptop.")
    return {"text": "\n".join(lines) + "\n", "info": info, "warnings": warnings, "name": f"{prefix}.scf.in"}
