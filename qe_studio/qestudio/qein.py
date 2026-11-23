"""Reading, validating and editing pw.x input files.

Edits are done on the text (not by re-serialising), so the user's own formatting and
comments survive in the generated inputs.
"""
from __future__ import annotations

import io
import math
import re
from collections import OrderedDict

import numpy as np

CARD_NAMES = (
    "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "ADDITIONAL_K_POINTS",
    "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES",
    "SOLVENTS", "HUBBARD",
)
CARD_RE = re.compile(r"^\s*(" + "|".join(CARD_NAMES) + r")\b(.*)$", re.I)
NL_ORDER = ["control", "system", "electrons", "ions", "cell", "fcp", "rism"]


# ----------------------------------------------------------------------------- text fixes
def normalize_text(text: str) -> tuple[str, list[str]]:
    """Fix Windows line endings and a missing final newline. Returns (text, notes)."""
    notes = []
    if "\r" in text:
        notes.append("Converted Windows (CRLF) line endings to Unix (LF).")
        text = text.replace("\r\n", "\n").replace("\r", "\n")
    if text.startswith("﻿"):
        text = text[1:]
        notes.append("Removed UTF-8 byte-order mark.")
    if not text.endswith("\n"):
        text += "\n"
        notes.append("Added the missing newline after the last line (QE fails to read the last namelist otherwise).")
    return text, notes


def _strip_comment(line: str) -> str:
    out, quote = [], None
    for ch in line:
        if quote:
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
        elif ch in "!#":
            break
        out.append(ch)
    return "".join(out)


def _split_commas(s: str) -> list[str]:
    parts, cur, quote = [], [], None
    for ch in s:
        if quote:
            if ch == quote:
                quote = None
            cur.append(ch)
        elif ch in "'\"":
            quote = ch
            cur.append(ch)
        elif ch == ",":
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return [p for p in parts if p.strip()]


def fortran_value(raw: str):
    v = raw.strip()
    if len(v) >= 2 and v[0] in "'\"" and v[-1] == v[0]:
        return v[1:-1]
    low = v.lower()
    if low in (".true.", "t", ".t.", "true"):
        return True
    if low in (".false.", "f", ".f.", "false"):
        return False
    try:
        return int(v)
    except ValueError:
        pass
    try:
        return float(low.replace("d", "e"))
    except ValueError:
        return v


def parse_input(text: str) -> dict:
    """Parse namelists and cards. Keys are lower-case; indexed keys keep their index: celldm(1)."""
    namelists: "OrderedDict[str, OrderedDict]" = OrderedDict()
    cards: "OrderedDict[str, dict]" = OrderedDict()
    cur_nl = None
    cur_card = None

    def eat(nl, s):
        for item in _split_commas(s):
            if "=" in item:
                k, v = item.split("=", 1)
                namelists[nl][re.sub(r"\s+", "", k).lower()] = fortran_value(v)

    for raw in text.split("\n"):
        line = _strip_comment(raw).strip()
        if cur_nl is not None:
            if not line:
                continue
            ended = False
            if line.endswith("/") and line.count("'") % 2 == 0 and line.count('"') % 2 == 0:
                line, ended = line[:-1].strip(), True
            eat(cur_nl, line)
            if ended:
                cur_nl = None
            continue
        m = re.match(r"^&(\w+)", line)
        if m:
            cur_nl = m.group(1).lower()
            namelists[cur_nl] = OrderedDict()
            rest = line[m.end():].strip()
            if rest.endswith("/"):
                eat(cur_nl, rest[:-1])
                cur_nl = None
            else:
                eat(cur_nl, rest)
            continue
        m = CARD_RE.match(line)
        if m:
            name = m.group(1).upper()
            opt = m.group(2).strip().strip("{}()").strip().lower()
            cur_card = {"option": opt, "lines": []}
            cards[name] = cur_card
            continue
        if cur_card is not None and line:
            cur_card["lines"].append(line)
    return {"namelists": namelists, "cards": cards}


def get(parsed: dict, nl: str, key: str, default=None):
    return parsed["namelists"].get(nl, {}).get(key.lower(), default)


# ----------------------------------------------------------------------------- text editing
def _fmt(value) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        return f"'{value}'"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def _namelist_span(text: str, nl: str):
    m = re.search(r"(?im)^[ \t]*&" + re.escape(nl) + r"\b.*$", text)
    if not m:
        return None
    end = re.compile(r"(?m)^[ \t]*/[ \t]*$").search(text, m.end())
    if not end:
        return None
    return m.start(), m.end(), end.start()


def set_param(text: str, nl: str, key: str, value) -> str:
    """Set key=value in namelist nl (adds the namelist if it is missing)."""
    val = _fmt(value)
    span = _namelist_span(text, nl)
    if span is None:
        block = f"&{nl.upper()}\n  {key} = {val}\n/\n"
        # insert after the last namelist that precedes it in the standard order
        idx = NL_ORDER.index(nl.lower()) if nl.lower() in NL_ORDER else len(NL_ORDER)
        insert_at = 0
        for prev in NL_ORDER[:idx]:
            sp = _namelist_span(text, prev)
            if sp:
                endline = text.find("\n", sp[2])
                insert_at = len(text) if endline < 0 else endline + 1
        return text[:insert_at] + block + text[insert_at:]
    start, head_end, close = span
    body = text[head_end:close]
    pat = re.compile(r"(?i)(?<![\w(])" + re.escape(key) + r"\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s/!]+)")
    if pat.search(body):
        body = pat.sub(lambda m: f"{key} = {val}", body, count=1)
        return text[:head_end] + body + text[close:]
    return text[:head_end] + f"\n  {key} = {val}" + text[head_end:]


def remove_param(text: str, nl: str, key: str) -> str:
    span = _namelist_span(text, nl)
    if span is None:
        return text
    start, head_end, close = span
    body = text[head_end:close]
    pat = re.compile(r"(?i)(?<![\w(])" + re.escape(key) + r"\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s/!]+)\s*,?")
    body = pat.sub("", body)
    body = re.sub(r"\n[ \t]*\n", "\n", body)
    return text[:head_end] + body + text[close:]


def replace_card(text: str, card: str, new_block: str) -> str:
    """Replace a card (header + its lines) with new_block; appends if missing."""
    lines = text.split("\n")
    start = None
    for i, ln in enumerate(lines):
        m = CARD_RE.match(_strip_comment(ln))
        if m and m.group(1).upper() == card.upper():
            start = i
            break
    if not new_block.endswith("\n"):
        new_block += "\n"
    if start is None:
        return text.rstrip("\n") + "\n" + new_block
    end = len(lines)
    for j in range(start + 1, len(lines)):
        s = _strip_comment(lines[j]).strip()
        if CARD_RE.match(s) or s.startswith("&"):
            end = j
            break
    new_lines = lines[:start] + new_block.rstrip("\n").split("\n") + lines[end:]
    out = "\n".join(new_lines)
    return out if out.endswith("\n") else out + "\n"


# ----------------------------------------------------------------------------- structure
def species_list(parsed: dict) -> list[dict]:
    out = []
    card = parsed["cards"].get("ATOMIC_SPECIES")
    if not card:
        return out
    for ln in card["lines"]:
        p = ln.split()
        if len(p) >= 3:
            out.append({"label": p[0], "mass": fortran_value(p[1]), "pseudo": p[2]})
    return out


BOHR = 0.529177210903  # Angstrom


def _celldm(parsed) -> dict:
    """celldm(1..6) from either celldm(i) or A,B,C,cosAB,cosAC,cosBC (QE conventions)."""
    s = parsed["namelists"].get("system", {})
    cd = {i: float(s[f"celldm({i})"]) for i in range(1, 7) if f"celldm({i})" in s}
    if "a" in s:
        A = float(s["a"])
        cd[1] = A / BOHR
        if "b" in s:
            cd[2] = float(s["b"]) / A
        if "c" in s:
            cd[3] = float(s["c"]) / A
        ibrav = int(s.get("ibrav", 0))
        cosab, cosac, cosbc = (float(s.get(k, 0.0)) for k in ("cosab", "cosac", "cosbc"))
        if ibrav == 14:
            cd[4], cd[5], cd[6] = cosbc, cosac, cosab
        elif ibrav in (-12, -13):
            cd[5] = cosac
        elif ibrav in (5, -5, 12, 13):
            cd[4] = cosab if "cosab" in s else cosbc
    return cd


def lattice_vectors(ibrav: int, cd: dict) -> np.ndarray:
    """Port of QE latgen: rows = a1,a2,a3 in units of alat (celldm(1))."""
    b = cd.get(2, 1.0)
    c = cd.get(3, 1.0)
    r3 = math.sqrt(3.0)
    if ibrav == 1:
        v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    elif ibrav == 2:
        v = [[-.5, 0, .5], [0, .5, .5], [-.5, .5, 0]]
    elif ibrav == 3:
        v = [[.5, .5, .5], [-.5, .5, .5], [-.5, -.5, .5]]
    elif ibrav == -3:
        v = [[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]]
    elif ibrav == 4:
        v = [[1, 0, 0], [-.5, r3 / 2, 0], [0, 0, c]]
    elif ibrav in (5, -5):
        cg = cd[4]
        tx, ty, tz = math.sqrt((1 - cg) / 2), math.sqrt((1 - cg) / 6), math.sqrt((1 + 2 * cg) / 3)
        if ibrav == 5:
            v = [[tx, -ty, tz], [0, 2 * ty, tz], [-tx, -ty, tz]]
        else:
            ap = 1 / r3
            u, w = tz - 2 * math.sqrt(2) * ty, tz + math.sqrt(2) * ty
            v = [[ap * u, ap * w, ap * w], [ap * w, ap * u, ap * w], [ap * w, ap * w, ap * u]]
    elif ibrav == 6:
        v = [[1, 0, 0], [0, 1, 0], [0, 0, c]]
    elif ibrav == 7:
        v = [[.5, -.5, c / 2], [.5, .5, c / 2], [-.5, -.5, c / 2]]
    elif ibrav == 8:
        v = [[1, 0, 0], [0, b, 0], [0, 0, c]]
    elif ibrav == 9:
        v = [[.5, b / 2, 0], [-.5, b / 2, 0], [0, 0, c]]
    elif ibrav == -9:
        v = [[.5, -b / 2, 0], [.5, b / 2, 0], [0, 0, c]]
    elif ibrav == 91:
        v = [[1, 0, 0], [0, b / 2, -c / 2], [0, b / 2, c / 2]]
    elif ibrav == 10:
        v = [[.5, 0, c / 2], [.5, b / 2, 0], [0, b / 2, c / 2]]
    elif ibrav == 11:
        v = [[.5, b / 2, c / 2], [-.5, b / 2, c / 2], [-.5, -b / 2, c / 2]]
    elif ibrav == 12:
        g = cd[4]
        v = [[1, 0, 0], [b * g, b * math.sqrt(1 - g * g), 0], [0, 0, c]]
    elif ibrav == -12:
        be = cd[5]
        v = [[1, 0, 0], [0, b, 0], [c * be, 0, c * math.sqrt(1 - be * be)]]
    elif ibrav == 13:
        g = cd[4]
        v = [[.5, 0, -c / 2], [b * g, b * math.sqrt(1 - g * g), 0], [.5, 0, c / 2]]
    elif ibrav == -13:
        be = cd[5]
        v = [[.5, b / 2, 0], [-.5, b / 2, 0], [c * be, 0, c * math.sqrt(1 - be * be)]]
    elif ibrav == 14:
        ca, cb, cg = cd[4], cd[5], cd[6]
        sg = math.sqrt(1 - cg * cg)
        v = [[1, 0, 0], [b * cg, b * sg, 0],
             [c * cb, c * (ca - cb * cg) / sg,
              c * math.sqrt(1 + 2 * ca * cb * cg - ca * ca - cb * cb - cg * cg) / sg]]
    else:
        raise ValueError(f"ibrav = {ibrav} is not supported")
    return np.array(v, dtype=float)


def _floats(line: str, n: int) -> list[float]:
    return [float(x.lower().replace("d", "e")) for x in line.split()[:n]]


def read_structure(text: str):
    """Build an ase.Atoms (Angstrom) from a pw.x input, for any ibrav."""
    from ase import Atoms
    parsed = parse_input(text) if isinstance(text, str) else text
    s = parsed["namelists"].get("system", {})
    ibrav = int(s.get("ibrav", 0))
    cd = _celldm(parsed)
    if ibrav == 0:
        card = parsed["cards"].get("CELL_PARAMETERS")
        if not card:
            raise ValueError("ibrav = 0 but CELL_PARAMETERS is missing")
        vec = np.array([_floats(l, 3) for l in card["lines"][:3]])
        opt = card["option"] or ("alat" if 1 in cd else "bohr")
        if opt == "angstrom":
            cell = vec
        elif opt == "bohr":
            cell = vec * BOHR
        else:  # alat
            if 1 not in cd:
                raise ValueError("CELL_PARAMETERS alat needs celldm(1) or A")
            cell = vec * cd[1] * BOHR
        alat_A = cd.get(1, np.linalg.norm(cell[0]) / BOHR) * BOHR if opt == "alat" else np.linalg.norm(cell[0])
    else:
        if 1 not in cd:
            raise ValueError("celldm(1) (or A) is required when ibrav != 0")
        alat_A = cd[1] * BOHR
        cell = lattice_vectors(ibrav, cd) * alat_A
    pos_card = parsed["cards"].get("ATOMIC_POSITIONS")
    if not pos_card:
        raise ValueError("ATOMIC_POSITIONS card is missing")
    opt = pos_card["option"] or "alat"
    nat = int(s.get("nat", len(pos_card["lines"])))
    labels, xyz = [], []
    for ln in pos_card["lines"][:nat]:
        p = ln.split()
        labels.append(p[0])
        xyz.append([float(x.lower().replace("d", "e")) for x in p[1:4]])
    xyz = np.array(xyz)
    if opt == "crystal":
        pos = xyz @ cell
    elif opt == "bohr":
        pos = xyz * BOHR
    elif opt == "angstrom":
        pos = xyz
    elif opt == "alat":
        pos = xyz * alat_A
    else:
        raise ValueError(f"ATOMIC_POSITIONS {opt} is not supported by the viewer")
    from ase.data import chemical_symbols
    symbols = []
    for lb in labels:
        m = re.match(r"([A-Z][a-z]?)", lb.capitalize())
        sym = m.group(1) if m and m.group(1) in chemical_symbols else "X"
        symbols.append(sym)
    atoms = Atoms(symbols=symbols, positions=pos, cell=cell, pbc=True)
    atoms.info["labels"] = labels
    atoms.info["alat_A"] = float(alat_A)
    atoms.info["ibrav"] = ibrav
    return atoms


def kpoints_info(parsed: dict) -> dict:
    card = parsed["cards"].get("K_POINTS")
    if not card:
        return {"type": "tpiba", "text": "K_POINTS missing (QE default: Gamma point)"}
    opt = card["option"] or "tpiba"
    info = {"type": opt}
    if opt == "automatic" and card["lines"]:
        nums = [int(x) for x in card["lines"][0].split()[:6]]
        info["grid"] = nums[:3]
        info["shift"] = nums[3:6] if len(nums) >= 6 else [0, 0, 0]
        info["text"] = "{} x {} x {} (shift {} {} {})".format(*info["grid"], *info["shift"])
    elif opt == "gamma":
        info["text"] = "Gamma only"
    else:
        info["text"] = f"{opt}: {card['lines'][0] if card['lines'] else '?'} points listed"
    return info


# ----------------------------------------------------------------------------- k-path
GREEK = {"G": "Γ", "GAMMA": "Γ", "S0": "S₀", "Y1": "Y₁", "X1": "X₁", "L1": "L₁", "Z1": "Z₁",
         "M1": "M₁", "N1": "N₁", "H1": "H₁", "F1": "F₁", "P1": "P₁", "A1": "A₁", "C1": "C₁",
         "D1": "D₁", "R1": "R₁"}


def pretty_label(lbl: str) -> str:
    return GREEK.get(lbl.upper(), lbl)


def make_kpath(atoms, path: str | None = None, density: float = 40.0) -> dict:
    """Build a crystal_b K_POINTS card along a high-symmetry path.

    path: ASE notation, e.g. "LGXUG" or "GXWKGLUWLK,UX"; None = ASE default for the lattice.
    density: points per (2pi/alat) unit... measured per 1/Angstrom of path length.
    Returns dict with card text, special-point labels, their indices in the k list and jumps.
    """
    cell = atoms.cell
    lat = cell.get_bravais_lattice()
    bp = cell.bandpath(npoints=0) if not path else cell.bandpath(path=path, npoints=0)
    special = bp.special_points  # scaled (crystal) coords for *this* cell
    segments = [seg for seg in re.findall(r"[A-Z][a-z0-9]*|,", bp.path)]
    # group into continuous segments
    groups, cur = [], []
    for tok in segments:
        if tok == ",":
            if cur:
                groups.append(cur)
            cur = []
        else:
            cur.append(tok)
    if cur:
        groups.append(cur)
    rec = cell.reciprocal() * 2 * math.pi  # rows b_i in 1/Angstrom
    pts, weights, labels = [], [], []
    for gi, grp in enumerate(groups):
        for j, name in enumerate(grp):
            k = np.array(special[name], dtype=float)
            if j + 1 < len(grp):
                k2 = np.array(special[grp[j + 1]], dtype=float)
                L = np.linalg.norm((k2 - k) @ rec)
                w = max(2, int(round(L * density)))
            else:
                w = 1  # jump to the next segment (or end of path)
            pts.append(k)
            weights.append(w)
            labels.append(name)
    lines = [f"K_POINTS crystal_b", str(len(pts))]
    for k, w, name in zip(pts, weights, labels):
        lines.append(f"  {k[0]: .8f} {k[1]: .8f} {k[2]: .8f} {w:4d} ! {name}")
    idx = [int(sum(weights[:i])) for i in range(len(pts))]
    jumps = []  # pairs of consecutive special points with no path between them
    pos = 0
    for gi, grp in enumerate(groups[:-1]):
        pos += len(grp)
        jumps.append(idx[pos - 1])
    return {
        "card": "\n".join(lines) + "\n",
        "labels": labels,
        "indices": idx,
        "jumps": jumps,
        "nks": int(sum(weights[:-1]) + 1),
        "path": bp.path,
        "lattice": f"{lat.longname} ({lat.name})",
        "special_points": {k: list(map(float, v)) for k, v in special.items()},
    }


def kgrid_for_density(atoms, n_max: int) -> list[int]:
    """Monkhorst-Pack grid with n_max along the longest reciprocal vector, scaled for others."""
    b = np.linalg.norm(atoms.cell.reciprocal(), axis=1)
    return [max(1, int(round(n_max * bi / b.max()))) for bi in b]


# ----------------------------------------------------------------------------- inspection
def inspect(text: str) -> dict:
    """Parse + validate an scf input. Never raises for bad input: errors go in 'errors'."""
    fixed, notes = normalize_text(text)
    res = {"notes": notes, "warnings": [], "errors": []}
    try:
        parsed = parse_input(fixed)
    except Exception as e:  # pragma: no cover
        res["errors"].append(f"Could not parse the input: {e}")
        return res
    nls = parsed["namelists"]
    res["namelists"] = {k: dict(v) for k, v in nls.items()}
    if "control" not in nls and "system" not in nls:
        if "bands" in nls or "dos" in nls or "projwfc" in nls:
            res["errors"].append("This is a post-processing input (bands.x/dos.x/projwfc.x), not a pw.x input. Upload your scf input.")
        else:
            res["errors"].append("No &CONTROL / &SYSTEM namelist found: this does not look like a pw.x input.")
        return res
    calc = str(get(parsed, "control", "calculation", "scf")).lower()
    res["calculation"] = calc
    if calc not in ("scf", "relax", "vc-relax"):
        res["warnings"].append(f"calculation = '{calc}'. The pipeline will run it as 'scf' first.")
    res["prefix"] = get(parsed, "control", "prefix", "pwscf")
    sp = species_list(parsed)
    res["species"] = sp
    if not sp:
        res["errors"].append("ATOMIC_SPECIES card is missing.")
    if "ATOMIC_POSITIONS" not in parsed["cards"]:
        res["errors"].append("ATOMIC_POSITIONS card is missing.")
    ecutwfc = get(parsed, "system", "ecutwfc")
    ecutrho = get(parsed, "system", "ecutrho")
    res["ecutwfc"] = ecutwfc
    res["ecutrho"] = ecutrho if ecutrho is not None else (4 * ecutwfc if isinstance(ecutwfc, (int, float)) else None)
    if ecutwfc is None:
        res["errors"].append("ecutwfc is not set in &SYSTEM.")
    elif ecutrho is not None and ecutrho < 4 * ecutwfc:
        res["errors"].append(f"ecutrho ({ecutrho}) must be at least 4 x ecutwfc ({4*ecutwfc}).")
    res["kpoints"] = kpoints_info(parsed)
    occ = str(get(parsed, "system", "occupations", "fixed")).lower()
    res["occupations"] = occ
    res["nspin"] = int(get(parsed, "system", "nspin", 1) or 1)
    res["noncolin"] = bool(get(parsed, "system", "noncolin", False))
    res["nbnd"] = get(parsed, "system", "nbnd")
    for key in ("outdir", "pseudo_dir"):
        v = get(parsed, "control", key)
        res[key] = v
    res["notes"].append("Each run gets its own folder: pseudo_dir is set to its ./pseudos/ and outdir to a fast Linux scratch folder. Your original file is kept unchanged.")
    try:
        atoms = read_structure(fixed)
        cell = atoms.cell
        res["structure"] = {
            "formula": atoms.get_chemical_formula(mode="hill"),
            "natoms": len(atoms),
            "symbols": atoms.get_chemical_symbols(),
            "positions": atoms.get_positions().round(6).tolist(),
            "cell": cell.array.round(6).tolist(),
            "abc": [round(x, 5) for x in cell.cellpar()],
            "volume": round(atoms.get_volume(), 4),
        }
        try:
            lat = cell.get_bravais_lattice()
            res["structure"]["lattice"] = f"{lat.longname} ({lat.name})"
        except Exception:
            res["structure"]["lattice"] = "unknown"
        try:
            kp = make_kpath(atoms)
            res["kpath"] = {"path": kp["path"], "nks": kp["nks"], "lattice": kp["lattice"],
                            "points": sorted(kp["special_points"].keys())}
        except Exception as e:
            res["warnings"].append(f"Could not determine a high-symmetry k-path automatically ({e}). Enter one by hand.")
        kp = res["kpoints"]
        if kp.get("grid"):
            # ignore directions sampled by a single k-point (vacuum in slabs / molecules)
            dens = [kp["grid"][i] / (np.linalg.norm(cell.reciprocal()[i]) * 2 * math.pi)
                    for i in range(3) if kp["grid"][i] > 1] or [99.0]
            res["kpoints"]["spacing"] = round(1.0 / min(dens), 4)  # 1/Angstrom per division
            if min(dens) < 2.5 and occ == "smearing":
                res["warnings"].append("The k-point grid looks coarse for a metal; check k-point convergence.")
            res["nscf_grid"] = [min(2 * g, 24) for g in kp["grid"]]
        else:
            res["nscf_grid"] = kgrid_for_density(atoms, 12)
    except Exception as e:
        res["errors"].append(f"Could not build the structure from the input: {e}")
    if occ == "fixed" and res["nspin"] == 2:
        res["warnings"].append("nspin=2 with fixed occupations needs tot_magnetization; consider smearing.")
    return res
