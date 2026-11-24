#!/usr/bin/env bash
# One-time setup for QE Studio inside WSL Ubuntu.
# Installs Miniforge (if missing) and a conda env "qe" with Quantum ESPRESSO 7.x from conda-forge.
# Never uses apt's quantum-espresso (QE 6.7 on Ubuntu 24.04 crashes in md5_from_file).
set -euo pipefail

MF="$HOME/miniforge3"
if [ ! -x "$MF/bin/conda" ]; then
  if ! command -v curl >/dev/null; then
    echo "ERROR: curl is missing. Install it in Ubuntu, then start again:"
    echo "       sudo apt update && sudo apt install -y curl"
    exit 1
  fi
  echo ">>> Installing Miniforge into $MF"
  cd "$HOME"
  curl -fL -o Miniforge3-Linux-x86_64.sh \
    https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash Miniforge3-Linux-x86_64.sh -b -p "$MF"
  "$MF/bin/conda" init bash >/dev/null
fi

if ! "$MF/bin/conda" env list | grep -q '^qe '; then
  echo ">>> Creating conda env 'qe' (Quantum ESPRESSO + Python tools)"
  "$MF/bin/conda" create -n qe -c conda-forge -y qe python numpy scipy matplotlib ase flask plotly
else
  echo ">>> Env 'qe' exists; making sure QE, Python and the packages are present"
  "$MF/bin/conda" install -n qe -c conda-forge -y qe python numpy scipy matplotlib ase flask plotly
fi

if dpkg -l 2>/dev/null | grep -q '^ii  quantum-espresso'; then
  echo "WARNING: apt package quantum-espresso is installed. It is broken on Ubuntu 24.04."
  echo "         Remove it with: sudo apt remove quantum-espresso"
fi

BIN="$MF/envs/qe/bin"
echo ">>> Installed tools:"
ls "$BIN"/pw.x "$BIN"/bands.x "$BIN"/dos.x "$BIN"/projwfc.x "$BIN"/mpirun
echo "" | OMP_NUM_THREADS=1 "$BIN/pw.x" 2>/dev/null | grep -m1 "Program PWSCF" || true
echo ">>> Setup complete. Start the app with: bash start.sh"
