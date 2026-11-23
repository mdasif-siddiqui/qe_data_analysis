#!/usr/bin/env bash
# Start QE Studio inside WSL. Open http://localhost:8050 in your Windows browser.
cd "$(dirname "$0")"
ENV="$HOME/miniforge3/envs/qe"
if [ ! -x "$ENV/bin/pw.x" ]; then
  echo "Quantum ESPRESSO env not found. Running setup first..."
  bash setup_wsl.sh || exit 1
fi
export PATH="$ENV/bin:$PATH"
export OMP_NUM_THREADS=1
exec python app.py
