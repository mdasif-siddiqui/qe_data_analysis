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
# PSlibrary valence tags (which semicore states are included; "l" = nonlinear core correction). Elements use
# different ones: Li/Be "s"/"sl", Mg "spnl", most p-block "n", transition metals "spn"/"dn" ... "" = H, He.
PSL_TAGS = ["n", "dn", "spn", "spdn", "s", "sl", "nl", "dnl", "spnl", "spdnl", "sp", "spl", ""]


def _psl(kind: str) -> list[str]:
    return [f"{{el}}.pbe-{t + '-' if t else ''}{kind}_psl.1.0.0.UPF" for t in PSL_TAGS]


CANDIDATES = {
    "pbe": _psl("rrkjus") + _psl("kjpaw") + ["{el}.pbe-rrkjus.UPF"],
    "pbe-paw": _psl("kjpaw"),
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


def _formula_counts(text: str) -> dict:
    return {el: float(n) if n else 1.0 for el, n in re.findall(r"([A-Z][a-z]?)\s*([\d.]*)", text)}


def _origin_setting(cif_text: str) -> tuple[str, str | None]:
    """ASE only honours origin choice 2 / rhombohedral axes via _symmetry_space_group_setting, but most CIFs
    write it as a suffix of the symbol ('F d -3 m :2', 'R -3 m :R'). Without symmetry operations in the file,
    ASE then uses setting 1 and silently builds a wrong structure (diamond with C at 1/8: 16 atoms, not 8).
    Returns the (possibly amended) CIF text and a note."""
    if re.search(r"_(symmetry_equiv_pos|space_group_symop)_\w*xyz", cif_text, re.I):
        return cif_text, None                      # explicit operations: the setting is in the operations
    if re.search(r"_symmetry_space_group_setting", cif_text, re.I):
        return cif_text, None
    m = re.search(r"^\s*_(?:symmetry_space_group_name_H-M|space_group_name_H-M_alt)\s+['\"]?([^'\"\n]+)", cif_text,
                  re.I | re.M)
    if not m:
        return cif_text, None
    suf = re.search(r":\s*([12RH])\s*$", m.group(1).strip(), re.I)
    if not suf:
        return cif_text, None
    setting = {"1": 1, "2": 2, "H": 1, "R": 2}[suf.group(1).upper()]
    if setting == 1:
        return cif_text, None
    what = "rhombohedral axes" if suf.group(1).upper() == "R" else "origin choice 2"
    eol = cif_text.find("\n", m.end())
    eol = len(cif_text) if eol < 0 else eol
    return (cif_text[:eol] + f"\n_symmetry_space_group_setting {setting}" + cif_text[eol:],
            f"The space group symbol '{m.group(1).strip()}' asks for {what}; it was applied.")


def _cif_checks(cif_text: str, atoms) -> list[str]:
    """Sanity checks on the structure built from the CIF (catches wrong symmetry/origin, overlaps, disorder)."""
    out: list[str] = []
    try:
        from ase.io.cif import parse_cif
        block = next(iter(parse_cif(io.StringIO(cif_text))))
    except Exception:
        block = None
    if block is not None:
        occ = block.get("_atom_site_occupancy")
        if occ is not None:
            vals = []
            for v in occ:
                try:
                    vals.append(float(str(v).split("(")[0]))
                except ValueError:
                    pass
            if any(v < 0.999 for v in vals):
                out.append("The CIF lists partial site occupancies (disorder). Each site is treated as fully occupied; "
                           "an ordered supercell is needed for a realistic calculation.")
        fsum = block.get("_chemical_formula_sum")
        if fsum:
            want = _formula_counts(str(fsum))
            have = _formula_counts(atoms.get_chemical_formula(mode="hill"))
            if want and set(want) == set(have):
                ref = next(iter(want))
                rw = {k: v / want[ref] for k, v in want.items()}
                rh = {k: v / have[ref] for k, v in have.items()}
                if any(abs(rw[k] - rh[k]) > 0.02 * max(1.0, rw[k]) for k in want):
                    out.append(f"The structure built from the CIF has composition {atoms.get_chemical_formula(mode='hill')}, "
                               f"but the CIF says {fsum}. Sites were doubled or lost: the coordinates probably do not match "
                               "the symmetry operations (origin choice). Check the CIF before running.")
            elif want:
                out.append(f"The CIF formula ({fsum}) has different elements than the structure built from it.")
            # atom count from Z (catches wrong settings also for one-element or ratio-preserving errors)
            z = block.get("_cell_formula_units_z")
            try:
                z = float(str(z).split("(")[0]) if z is not None else None
            except ValueError:
                z = None
            if z and want and set(want) == set(have) and not out:
                expect = z * sum(want.values())
                if abs(expect - len(atoms)) > 0.5:
                    out.append(f"The CIF says Z = {z:g} formula units ({expect:g} atoms in the cell), but {len(atoms)} atoms "
                               "were built: the coordinates probably do not match the symmetry setting (origin choice). "
                               "Check the CIF before running.")
    if len(atoms) > 1:
        d = atoms.get_all_distances(mic=True)
        np.fill_diagonal(d, np.inf)
        if d.min() < 0.7:
            out.append(f"Two atoms are only {d.min():.2f} Å apart: overlapping sites, usually from duplicated positions "
                       "or a wrong symmetry setting.")
    return out


# Elements whose compounds are often magnetic (open 3d / 4f shells). "auto" spin polarisation starts these
# ferromagnetically; non-magnetic compounds of them (e.g. low-spin Fe2+) simply converge to zero moment.
MAGNETIC = {"Cr", "Mn", "Fe", "Co", "Ni", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"}


def cif_to_scf(cif_text: str, functional: str = "pbe", primitive: bool = True, occupations: str = "smearing",
               kspacing: float = 0.25, prefix: str | None = None, download: bool = True, spin: str = "auto") -> dict:
    """Returns {"text": scf input, "info": {...}, "warnings": [...]}."""
    from ase.io import read
    from ase import Atoms
    warnings: list[str] = []
    if functional not in CANDIDATES:
        raise ValueError(f"Unknown functional {functional!r}")
    cif_text, note = _origin_setting(cif_text)
    if note:
        warnings.append(note)
    try:
        atoms = read(io.StringIO(cif_text), format="cif")
    except Exception as e:
        raise ValueError(f"Could not read the CIF file: {e}")
    if isinstance(atoms, list):
        atoms = atoms[0]
    warnings += _cif_checks(cif_text, atoms)
    info: dict = {"formula_cif": atoms.get_chemical_formula(mode="hill"), "natoms_cif": len(atoms)}
    try:
        import spglib
        cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
        ds = spglib.get_symmetry_dataset(cell, symprec=1e-3)
        if ds is not None:
            info["spacegroup"] = f"{ds.international} (No. {ds.number})"
        if primitive:
            # standardized primitive cell with idealized positions: rounded CIF values such as 0.33333
            # would otherwise break the symmetry slightly (QE then finds fewer symmetry operations)
            lat, pos, nums = spglib.standardize_cell(cell, to_primitive=True, no_idealize=False, symprec=1e-3)
            if len(nums) <= len(atoms):
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
    # The old norm-conserving files (pz-vbc, ...) understate their needs for hard first-row elements:
    # AlN with N.pz-vbc ("suggested 60 Ry") still changes by 82 mRy/cell and 130 kbar from 60 to 80 Ry and
    # needs about 100 Ry.
    hard = sorted(el for el in elements if el in ("B", "C", "N", "O", "F"))
    if nc_only and hard and ecutwfc < 100:
        ecutwfc = 100.0
        ecw = ecw or [ecutwfc]          # the cutoff is now chosen, not a blind guess
        warnings.append(f"Norm-conserving pseudopotentials for {', '.join(hard)} need about 100 Ry; ecutwfc was raised to "
                        "100 Ry (the files suggest less, which is not converged). Run the cutoff convergence test to check.")
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
    mag = [el for el in elements if el in MAGNETIC]
    polarised = spin == "on" or (spin == "auto" and bool(mag))
    info["spin"] = "spin-polarised (ferromagnetic start)" if polarised else "non-magnetic"
    if polarised:
        lines.append("  nspin = 2")
        for i, el in enumerate(elements, 1):
            if el in MAGNETIC or spin == "on":
                lines.append(f"  starting_magnetization({i}) = {0.5 if el in MAGNETIC else 0.1}")
        warnings.append(f"Spin-polarised with a ferromagnetic start on {', '.join(mag) or 'all atoms'}. Antiferromagnetic or "
                        "other orders need a supercell with separate species for spin-up and spin-down sites.")
        if occupations != "smearing":
            lines.append("  tot_magnetization = 0")
            warnings.append("Fixed occupations with nspin = 2 need a fixed total moment: tot_magnetization = 0 was set "
                            "(non-magnetic result). For a magnet, choose 'metal / unknown (smearing)'.")
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
    if functional == "lda-nc":
        warnings.append("For IR/Raman (and phonons) choose 'relax first: vc-relax' under SCF: at the experimental "
                        "lattice LDA gives a large pressure and frequencies that are several percent off.")
    if len(atoms) > 40:
        warnings.append(f"{len(atoms)} atoms: this will be slow on a laptop.")
    return {"text": "\n".join(lines) + "\n", "info": info, "warnings": warnings, "name": f"{prefix}.scf.in"}
