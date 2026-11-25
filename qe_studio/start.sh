#!/usr/bin/env bash
# Start QE Studio inside WSL. Open http://localhost:8050 in your Windows browser.
cd "$(dirname "$0")" || exit 1

if ! grep -qiE 'ubuntu|debian' /etc/os-release 2>/dev/null; then
  echo "WARNING: this WSL system is not Ubuntu ($(. /etc/os-release 2>/dev/null; echo "${PRETTY_NAME:-unknown}"))."
  echo "         Make Ubuntu the default in PowerShell:  wsl --set-default Ubuntu"
fi

ENV="$HOME/miniforge3/envs/qe"
PY="$ENV/bin/python"

# Setup is needed if QE *or* Python and the app's packages are missing
# (an existing 'qe' env made with only 'conda create -n qe qe' has pw.x but no Python).
need_setup=0
[ -x "$ENV/bin/pw.x" ] || need_setup=1
{ [ -x "$PY" ] && "$PY" -c "import flask, plotly, ase, numpy, spglib" 2>/dev/null; } || need_setup=1

if [ "$need_setup" = 1 ]; then
  echo "Quantum ESPRESSO or its Python packages are missing. Running setup (needs internet, about 5-10 minutes)..."
  bash setup_wsl.sh || { echo; echo "Setup failed - see the messages above."; exit 1; }
fi

if [ ! -x "$PY" ]; then
  echo "ERROR: $PY is still missing after setup. Please send the messages above."
  exit 1
fi

export PATH="$ENV/bin:$PATH"
export OMP_NUM_THREADS=1
exec "$PY" app.py
