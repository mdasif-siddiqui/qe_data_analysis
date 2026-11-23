"""Settings from environment variables.

Local use needs none of these. For an internet-facing deployment set QE_STUDIO_PUBLIC=1:
that disables server-side folder access, Explorer, and deleting runs, and enforces the limits below.
"""
from __future__ import annotations

import os


def _int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, default))
    except ValueError:
        return default


PUBLIC = os.environ.get("QE_STUDIO_PUBLIC", "0").lower() in ("1", "true", "yes")
PASSWORD = os.environ.get("QE_STUDIO_PASSWORD", "")      # if set: HTTP basic auth on every page
HOST = os.environ.get("QE_STUDIO_HOST", "127.0.0.1")
PORT = _int("QE_STUDIO_PORT", 8050)
QE_BIN = os.environ.get("QE_BIN", "")                     # folder with pw.x; default: ~/miniforge3/envs/qe/bin or PATH

# limits (always available; defaults are generous locally, strict in public mode)
MAX_ATOMS = _int("QE_STUDIO_MAX_ATOMS", 40 if PUBLIC else 100000)
MAX_ECUTWFC = _int("QE_STUDIO_MAX_ECUTWFC", 100 if PUBLIC else 100000)
MAX_KPOINTS = _int("QE_STUDIO_MAX_KPOINTS", 4096 if PUBLIC else 10**9)      # n1*n2*n3 of any grid
MAX_STEP_SECONDS = _int("QE_STUDIO_MAX_STEP_SECONDS", 900 if PUBLIC else 0)  # 0 = no limit
MAX_QUEUE = _int("QE_STUDIO_MAX_QUEUE", 5 if PUBLIC else 1000)
MAX_UPLOAD_MB = _int("QE_STUDIO_MAX_UPLOAD_MB", 20 if PUBLIC else 500)
KEEP_DAYS = _int("QE_STUDIO_KEEP_DAYS", 3 if PUBLIC else 0)                 # delete older runs; 0 = keep
MAX_NP = _int("QE_STUDIO_MAX_NP", 0)                                         # 0 = all available cores
