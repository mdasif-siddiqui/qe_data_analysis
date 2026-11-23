"""Find, validate and download pseudopotentials."""
from __future__ import annotations

import os
import re
import shutil
import urllib.request
from pathlib import Path

LIBRARY = Path(__file__).resolve().parent.parent / "pseudo_library"
QE_URL = "https://pseudopotentials.quantum-espresso.org/upf_files/{name}"


NAME_RE = re.compile(r"^[A-Za-z0-9][\w.+\-]{0,120}$")


def safe_name(name: str) -> str:
    """Pseudopotential file name without any directory part; rejects odd names."""
    base = str(name)
    if not NAME_RE.match(base):  # no directory parts: QE would follow them relative to pseudo_dir
        raise ValueError(f"Invalid pseudopotential file name: {name!r}")
    return base


def is_valid_upf(path: Path) -> tuple[bool, str]:
    try:
        head = path.read_bytes()[:4000].decode("utf-8", "replace")
    except OSError as e:
        return False, str(e)
    low = head.lower()
    if "<html" in low or "<!doctype" in low:
        return False, "file is an HTML page (failed download), not a pseudopotential"
    if "<upf" in low or "<pp_info" in low or "<pp_header" in low:
        return True, "ok"
    return False, "file does not look like a UPF pseudopotential"


def z_valence(path: Path) -> float | None:
    txt = path.read_text(errors="replace")[:200000]
    m = re.search(r'z_valence\s*=\s*"\s*([-+0-9.eEdD]+)', txt)
    if m:
        return float(m.group(1).lower().replace("d", "e"))
    m = re.search(r"([-+0-9.eE]+)\s+Z valence", txt)
    if m:
        return float(m.group(1))
    return None


def search_dirs(extra: list[Path] | None = None) -> list[Path]:
    dirs = list(extra or [])
    dirs.append(LIBRARY)
    for env in ("ESPRESSO_PSEUDO",):
        if os.environ.get(env):
            dirs.append(Path(os.environ[env]))
    dirs.append(Path.home() / "espresso" / "pseudo")
    return [d for d in dirs if d.is_dir()]


def find(name: str, extra: list[Path] | None = None) -> Path | None:
    try:
        name = safe_name(name)
    except ValueError:
        return None
    for d in search_dirs(extra):
        p = d / name
        if p.is_file() and is_valid_upf(p)[0]:
            return p
    return None


def download(name: str) -> Path:
    name = safe_name(name)
    LIBRARY.mkdir(parents=True, exist_ok=True)
    dest = LIBRARY / name
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(QE_URL.format(name=name), headers={"User-Agent": "qe-studio"})
    with urllib.request.urlopen(req, timeout=60) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f)
    ok, why = is_valid_upf(tmp)
    if not ok:
        tmp.unlink(missing_ok=True)
        raise RuntimeError(f"{name}: download from the QE library failed ({why}). Upload the file yourself.")
    tmp.replace(dest)
    return dest


def ensure(names: list[str], dest_dir: Path, extra: list[Path] | None = None,
           allow_download: bool = True, log=print) -> dict:
    """Copy each pseudo into dest_dir. Returns {name: {"source":..., "z_valence":...}}."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    info = {}
    for name in names:
        name = safe_name(name)
        target = dest_dir / name
        if target.is_file() and is_valid_upf(target)[0]:
            src = "uploaded"
        else:
            p = find(name, extra)
            if p is None:
                if not allow_download:
                    raise RuntimeError(f"Pseudopotential {name} not found. Upload it or enable automatic download.")
                log(f"Downloading {name} from pseudopotentials.quantum-espresso.org ...")
                p = download(name)
                src = "downloaded from the QE library"
            else:
                src = str(p)
            shutil.copy2(p, target)
        # keep a copy in the library so later runs need no download/upload
        LIBRARY.mkdir(parents=True, exist_ok=True)
        if not (LIBRARY / name).exists():
            shutil.copy2(target, LIBRARY / name)
        info[name] = {"source": src, "z_valence": z_valence(target)}
    return info


def library_listing() -> list[str]:
    if not LIBRARY.is_dir():
        return []
    return sorted(p.name for p in LIBRARY.iterdir() if p.is_file() and not p.name.endswith(".part"))
