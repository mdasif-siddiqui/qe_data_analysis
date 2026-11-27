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


def _psl(kind: str, rel: bool = False) -> list[str]:
    return [f"{{el}}.{'rel-' if rel else ''}pbe-{t + '-' if t else ''}{kind}_psl.1.0.0.UPF" for t in PSL_TAGS]


CANDIDATES = {
    "pbe": _psl("rrkjus") + _psl("kjpaw") + ["{el}.pbe-rrkjus.UPF"],
    "pbe-paw": _psl("kjpaw"),
    "lda-nc": ["{el}.pz-vbc.UPF", "{el}.pz-bhs.UPF", "{el}.pz-mt.UPF", "{el}.pz-n-vbc.UPF", "{el}.pz-hgh.UPF"],
    # fully relativistic PSlibrary files: spin-orbit coupling (noncolin + lspinorb)
    "pbe-so": _psl("rrkjus", rel=True) + _psl("kjpaw", rel=True),
}
FUNCTIONALS = {"pbe": "PBE (ultrasoft/PAW, PSlibrary)", "pbe-paw": "PBE PAW (PSlibrary; GIPAW-ready for NMR)",
               "lda-nc": "LDA, norm-conserving (needed for Raman and UV-Vis)",
               "pbe-so": "PBE with spin-orbit coupling (fully relativistic PSlibrary)"}
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


def _site_labels(cif_text: str, atoms) -> list[str] | None:
    """CIF site label of every atom ASE built (ASE keeps the site index in the 'spacegroup_kinds' array)."""
    if "spacegroup_kinds" not in atoms.arrays:
        return None
    try:
        from ase.io.cif import parse_cif
        block = next(iter(parse_cif(io.StringIO(cif_text))))
        labels = [str(x) for x in block.get("_atom_site_label")]
    except Exception:
        return None
    kinds = atoms.arrays["spacegroup_kinds"]
    if len(labels) <= int(max(kinds)):
        return None
    return [labels[int(k)] for k in kinds]


# Elements whose compounds are often magnetic (open 3d / 4f shells). "auto" spin polarisation starts these
# ferromagnetically; non-magnetic compounds of them (e.g. low-spin Fe2+) simply converge to zero moment.
MAGNETIC = {"Cr", "Mn", "Fe", "Co", "Ni", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb"}

# Hubbard U (eV) of the Materials Project (GGA+U, fitted to oxide formation enthalpies), used by
# hubbard="auto" for oxides and fluorides. Any element can be given by hand instead ("Ni=6.2, Fe=5").
HUBBARD_U = {"V": 3.25, "Cr": 3.7, "Mn": 3.9, "Fe": 5.3, "Co": 3.32, "Ni": 6.2, "Mo": 4.38, "W": 6.2}


def _hubbard_shell(el: str) -> str | None:
    from ase.data import atomic_numbers
    z = atomic_numbers.get(el, 0)
    if 21 <= z <= 30:
        return "3d"
    if 39 <= z <= 48:
        return "4d"
    if 72 <= z <= 80:
        return "5d"
    if 57 <= z <= 71:
        return "4f"
    if 89 <= z <= 103:
        return "5f"
    return None


def _has_wfc(upf_path, shell: str) -> bool:
    """True if the UPF has an atomic wavefunction labelled like '3D' (needed as the Hubbard projector)."""
    try:
        txt = open(upf_path, errors="replace").read()
    except OSError:
        return True   # cannot check: let pw.x decide
    return re.search(r'label\s*=\s*"\s*' + shell.upper() + r'\s*"', txt, re.I) is not None or \
        re.search(r"^\s*" + shell.upper() + r"\s+\d", txt, re.M) is not None


def parse_hubbard(spec: str) -> dict[str, float]:
    """'Ni=6.2, Fe = 5' -> {'Ni': 6.2, 'Fe': 5.0}."""
    out = {}
    for el, v in re.findall(r"([A-Z][a-z]?)\s*[=:]\s*([0-9.]+)", spec or ""):
        out[el] = float(v)
    return out


def cif_to_scf(cif_text: str, functional: str = "pbe", primitive: bool = True, occupations: str = "smearing",
               kspacing: float = 0.25, prefix: str | None = None, download: bool = True, spin: str = "auto",
               hubbard: str = "auto") -> dict:
    """Returns {"text": scf input, "info": {...}, "warnings": [...]}.

    hubbard: "auto" (Materials Project U for V, Cr, Mn, Fe, Co, Ni, Mo, W in oxides and fluorides), "off",
    or explicit values such as "Ni=6.2, Fe=5.3" (eV)."""
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

    # species per atom: normally the element; for antiferromagnets each CIF site label of a magnetic element
    # becomes its own species (Ni1, Ni2, ...), with alternating starting moments
    symbols = atoms.get_chemical_symbols()
    species = list(symbols)
    afm_sign: dict[str, int] = {}
    if spin == "afm":
        site_of = _site_labels(cif_text, atoms)
        for el in dict.fromkeys(symbols):
            if el not in MAGNETIC or site_of is None:
                continue
            labs = list(dict.fromkeys(site_of[i] for i, s in enumerate(symbols) if s == el))
            if len(labs) < 2:
                continue
            for j, lab in enumerate(labs):
                afm_sign[f"{el}{j + 1}"] = 1 if j % 2 == 0 else -1
            for i, s in enumerate(symbols):
                if s == el:
                    species[i] = f"{el}{labs.index(site_of[i]) + 1}"
        if not afm_sign:
            warnings.append("Antiferromagnetic needs the magnetic cell in the CIF, with the spin-up and spin-down sites of "
                            "the magnetic element under different labels (e.g. Ni1 and Ni2). This CIF has one label per "
                            "magnetic element, so a ferromagnetic start was used instead.")
    kinds_list = list(dict.fromkeys(species))
    type_id = {sp: i + 1 for i, sp in enumerate(kinds_list)}
    sp_el = {sp: el for sp, el in zip(species, symbols)}
    try:
        import spglib
        cell = (atoms.cell.array, atoms.get_scaled_positions(), [type_id[s] for s in species])
        ds = spglib.get_symmetry_dataset(cell, symprec=1e-3)
        if ds is not None:
            info["spacegroup"] = f"{ds.international} (No. {ds.number})" + (" (magnetic cell)" if afm_sign else "")
        if primitive:
            # standardized primitive cell with idealized positions: rounded CIF values such as 0.33333
            # would otherwise break the symmetry slightly (QE then finds fewer symmetry operations).
            # Species ids (not atomic numbers) keep spin-up and spin-down sites apart.
            lat, pos, nums = spglib.standardize_cell(cell, to_primitive=True, no_idealize=False, symprec=1e-3)
            if len(nums) <= len(atoms):
                inv = {v: k for k, v in type_id.items()}
                species = [inv[int(n)] for n in nums]
                atoms = Atoms(symbols=[sp_el[s] for s in species], cell=lat, scaled_positions=pos, pbc=True)
    except ImportError:
        warnings.append("spglib is not installed: the cell was not reduced and the space group is unknown.")
    info["natoms"] = len(atoms)
    info["formula"] = atoms.get_chemical_formula(mode="hill")
    elements = list(dict.fromkeys(atoms.get_chemical_symbols()))
    species_list = list(dict.fromkeys(species))
    from ase.data import atomic_masses, atomic_numbers

    # pseudopotentials and cutoffs
    ps, ecw, ecr, kinds, ps_path = {}, [], [], set(), {}
    for el in elements:
        name = choose_pseudo(el, functional)
        ps[el] = name
        ps_path[el] = None
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
            ps_path[el] = path
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
             f"  nat = {len(atoms)}", f"  ntyp = {len(species_list)}", f"  ecutwfc = {ecutwfc:g}", f"  ecutrho = {ecutrho:g}"]
    if occupations == "smearing":
        lines += ["  occupations = 'smearing'", "  smearing = 'mv'", "  degauss = 0.01"]
    mag = [el for el in elements if el in MAGNETIC]
    if functional == "pbe-so":
        # spin-orbit: two-component spinors. Magnetic order then needs non-collinear starting moments (angle1/2),
        # which this converter does not set up: SOC inputs are made non-magnetic.
        lines += ["  noncolin = .true.", "  lspinorb = .true."]
        info["spin"] = "spin-orbit coupling (noncolin, lspinorb), non-magnetic"
        if mag and spin != "off":
            warnings.append(f"Spin-orbit coupling: {', '.join(mag)} is usually magnetic, but magnetic SOC inputs "
                            "(starting_magnetization with angle1/angle2) are not set up by the converter; this input is "
                            "non-magnetic. Add them by hand if needed.")
        spin = "off"
    polarised = spin in ("on", "afm") or (spin == "auto" and bool(mag))
    info["spin"] = info["spin"] if functional == "pbe-so" else ("spin-polarised (antiferromagnetic start: " + ", ".join(f"{s} {'+' if v > 0 else '−'}" for s, v in afm_sign.items()) + ")"
                    if afm_sign else "spin-polarised (ferromagnetic start)" if polarised else "non-magnetic")
    if polarised:
        lines.append("  nspin = 2")
        for i, sp in enumerate(species_list, 1):
            el = sp_el[sp]
            if sp in afm_sign:
                lines.append(f"  starting_magnetization({i}) = {0.5 * afm_sign[sp]:g}")
            elif el in MAGNETIC or spin == "on":
                lines.append(f"  starting_magnetization({i}) = {0.5 if el in MAGNETIC else 0.1}")
        if afm_sign:
            warnings.append("Antiferromagnetic start: " + ", ".join(f"{s} {'up' if v > 0 else 'down'}" for s, v in afm_sign.items())
                            + ". Check the moments in the results; a calculation can still end in another magnetic state.")
        else:
            warnings.append(f"Spin-polarised with a ferromagnetic start on {', '.join(mag) or 'all atoms'}. For an "
                            "antiferromagnet choose Spin: antiferromagnetic and give the magnetic cell with separate site labels.")
        if occupations != "smearing":
            lines.append("  tot_magnetization = 0")
            if not afm_sign:
                warnings.append("Fixed occupations with nspin = 2 need a fixed total moment: tot_magnetization = 0 was set "
                                "(non-magnetic result). For a ferromagnet, choose 'metal / unknown (smearing)'.")
    lines += ["/", "&ELECTRONS", "  conv_thr = 1.0d-9", "  mixing_beta = 0.5", "/", "ATOMIC_SPECIES"]
    for sp in species_list:
        el = sp_el[sp]
        mass = atomic_masses[atomic_numbers[el]]
        lines.append(f"  {sp:3s} {mass:9.4f} {ps[el] or el + '.UPF'}")
    lines.append("CELL_PARAMETERS angstrom")
    for v in atoms.cell.array:
        lines.append("  " + " ".join(f"{x:14.9f}" for x in v))
    lines.append("ATOMIC_POSITIONS crystal")
    for s, p in zip(species, atoms.get_scaled_positions(wrap=True)):
        lines.append(f"  {s:3s} " + " ".join(f"{x:13.9f}" for x in p))
    lines += ["K_POINTS automatic", "  {} {} {} 0 0 0".format(*grid)]

    # DFT+U (Dudarev, ortho-atomic projectors, QE 7 HUBBARD card)
    hub = str(hubbard or "auto").strip()
    if hub.lower() == "auto":
        anion = any(el in ("O", "F") for el in elements)
        u_of = {el: HUBBARD_U[el] for el in elements if el in HUBBARD_U} if anion else {}
        source = "Materials Project values"
    elif hub.lower() in ("off", "none", "0", ""):
        u_of, source = {}, ""
    else:
        u_of, source = parse_hubbard(hub), "given"
        if not u_of:
            warnings.append(f"Could not read the Hubbard U values {hub!r}; write them like 'Ni=6.2, Fe=5.3'. No +U was used.")
        for el in [e for e in u_of if e not in elements]:
            warnings.append(f"Hubbard U for {el} ignored: there is no {el} in this structure.")
            u_of.pop(el)
    hub_lines, used = [], []
    for sp in species_list:
        el = sp_el[sp]
        if el not in u_of or u_of[el] <= 0:
            continue
        shell = _hubbard_shell(el)
        if shell is None:
            warnings.append(f"No d or f shell known for {el}: Hubbard U not applied to it.")
            continue
        if ps_path.get(el) and not _has_wfc(ps_path[el], shell):
            warnings.append(f"The pseudopotential {ps[el]} has no {shell} atomic wavefunction, which DFT+U needs as "
                            f"projector: U on {el} was not applied.")
            continue
        hub_lines.append(f"  U {sp}-{shell} {u_of[el]:g}")
        used.append(f"{sp}-{shell} {u_of[el]:g} eV")
    if hub_lines:
        lines += ["HUBBARD (ortho-atomic)"] + hub_lines
        info["hubbard"] = "DFT+U (" + source + "): " + ", ".join(dict.fromkeys(used))
        warnings.append("DFT+U: " + ", ".join(dict.fromkeys(used)) + f" ({source}, ortho-atomic projectors). U shifts "
                        "band gaps, moments and energies; compare only calculations with the same U. Phonons, IR/Raman "
                        "and NMR with +U are not tested.")
    else:
        info["hubbard"] = "none"
    if occupations == "fixed":
        warnings.append("Fixed occupations: correct only for insulators and semiconductors. Use smearing for metals.")
    if functional == "lda-nc":
        warnings.append("For IR/Raman (and phonons) choose 'relax first: vc-relax' under SCF: at the experimental "
                        "lattice LDA gives a large pressure and frequencies that are several percent off.")
    if len(atoms) > 40:
        warnings.append(f"{len(atoms)} atoms: this will be slow on a laptop.")
    return {"text": "\n".join(lines) + "\n", "info": info, "warnings": warnings, "name": f"{prefix}.scf.in"}
