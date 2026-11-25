"""IR and Raman spectra from a Gamma-point ph.x run (epsil + lraman) and dynmat.x."""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

IR_KMMOL = 42.2561          # 1 (D/Angstrom)^2/amu in km/mol
F = r"[-+]?\d*\.\d+(?:[eEdD][-+]?\d+)?"


def ph_gamma_input(prefix: str, outdir: str, raman: bool) -> str:
    return (f"phonons at Gamma\n&INPUTPH\n  prefix = '{prefix}'\n  outdir = '{outdir}'\n  fildyn = '{prefix}.dynG'\n"
            f"  tr2_ph = 1.0d-16\n  epsil = .true.\n  lraman = {'.true.' if raman else '.false.'}\n/\n0.0 0.0 0.0\n")


def dynmat_input(prefix: str, q_dir: list[float] | None = None) -> str:
    """dynmat.x writes eigenvectors to 'dynmat.out' by default; give them their own names so they never
    overwrite the program's screen output (which holds the IR/Raman table)."""
    q = "" if q_dir is None else f"\n  q(1) = {q_dir[0]}, q(2) = {q_dir[1]}, q(3) = {q_dir[2]}"
    tag = "_lo" if q_dir is not None else ""
    return (f"&INPUT\n  fildyn = '{prefix}.dynG'\n  asr = 'crystal'{q}\n  filout = '{prefix}.modes{tag}'\n"
            f"  filmol = '{prefix}.mold{tag}'\n  filxsf = '{prefix}.axsf{tag}'\n/\n")


def raman_allowed(xc: str | None, pseudo_types: list[str], occupations: str) -> tuple[bool, str]:
    """ph.x computes Raman tensors only for insulators with LDA and norm-conserving pseudopotentials."""
    if occupations not in ("fixed", ""):
        return False, "Raman (and IR) need an insulator with occupations = 'fixed'."
    if xc and re.search(r"\bPBE\b|\bPBX\b|PBC|B88|GGA|LYP|PW91|SCAN|HSE|B3LYP", xc.upper()):
        return False, f"Raman tensors in ph.x need LDA; this calculation uses {xc}. Use LDA norm-conserving pseudopotentials."
    if any(t not in ("NC", "SL") for t in pseudo_types):
        return False, "Raman tensors in ph.x need norm-conserving pseudopotentials (these are ultrasoft/PAW)."
    return True, ""


def parse_ph_dielectric(text: str) -> dict:
    res: dict = {}
    blocks = re.findall(r"Dielectric constant in cartesian axis\s*\n\s*\n((?:\s*\(.*\)\s*\n){3})", text)
    if blocks:
        rows = [[float(x) for x in re.findall(F, l)] for l in blocks[-1].strip().splitlines()]
        res["eps_inf"] = rows
        res["eps_inf_mean"] = float(np.trace(np.array(rows)) / 3)
    m = list(re.finditer(r"Effective charges \(d Force / dE\) in cartesian axis without acoustic sum rule applied \(asr\)\s*\n(.*?)\n\s*Effective charges Sum", text, re.S))
    if m:
        zs = []
        for am in re.finditer(r"atom\s+(\d+)\s+(\w+)\s+Mean Z\*:\s*(" + F + r")\s*\n((?:\s*E[xyz]\s*\(.*\)\s*\n){3})", m[-1].group(1) + "\n"):
            tensor = [[float(x) for x in re.findall(F, l)] for l in am.group(4).strip().splitlines()]
            zs.append({"atom": int(am.group(1)), "element": am.group(2), "mean": float(am.group(3)), "tensor": tensor})
        res["born"] = zs
        s = re.search(r"Effective charges Sum: Mean:\s*(" + F + ")", text[m[-1].start():])
        res["born_sum"] = float(s.group(1)) if s else None
    return res


def parse_dynmat(text: str) -> list[dict]:
    """Mode table printed by dynmat.x: frequency, IR intensity, Raman activity, depolarisation."""
    block = text.split("# mode", 1)[1] if "# mode" in text else ""
    modes = []
    for line in block.splitlines()[1:]:
        nums = re.findall(F, line)
        if not re.match(r"\s*\d+\s", line) or len(nums) < 3:
            if modes:
                break
            continue
        n = int(line.split()[0])
        vals = [float(x) for x in nums]
        mode = {"mode": n, "freq_cm1": vals[0], "ir": vals[2] * IR_KMMOL if len(vals) > 2 else None}
        if len(vals) >= 5:
            mode["raman"], mode["depol"] = vals[3], vals[4]
        modes.append(mode)
    return modes


def analyse(d: Path, prefix: str) -> dict | None:
    sd = d / "spectra"
    if not (sd / "ph.out").exists():
        return None
    ph = (sd / "ph.out").read_text(errors="replace")
    res = parse_ph_dielectric(ph)
    res["raman_computed"] = "Raman tensor" in ph
    if (sd / "dynmat.out").exists():
        res["modes"] = parse_dynmat((sd / "dynmat.out").read_text(errors="replace"))
    if (sd / "dynmat_lo.out").exists() and res.get("modes"):
        lo = parse_dynmat((sd / "dynmat_lo.out").read_text(errors="replace"))
        for m, l in zip(res["modes"], lo):
            m["freq_lo_cm1"] = l["freq_cm1"]
    modes = [m for m in res.get("modes", []) if m["freq_cm1"] > 10]
    if modes:
        irmax = max((m["ir"] or 0) for m in modes) or 1
        rmax = max((m.get("raman") or 0) for m in modes) or 1
        for m in modes:
            m["ir_active"] = (m["ir"] or 0) > 1e-3 * irmax and (m["ir"] or 0) > 1e-4
            m["raman_active"] = (m.get("raman") or 0) > 1e-3 * rmax and (m.get("raman") or 0) > 1e-4
        lo_split = [m["freq_lo_cm1"] - m["freq_cm1"] for m in modes if m.get("freq_lo_cm1") is not None]
        res["max_lo_to_cm1"] = max(lo_split) if lo_split else None
    return res
