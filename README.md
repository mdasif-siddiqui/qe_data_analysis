# QE Data Analysis

**QE Studio** is a local web app for running and analysing [Quantum ESPRESSO](https://www.quantum-espresso.org/)
calculations on Windows via WSL Ubuntu.

Give it a pw.x `scf.in` and it validates the input, runs SCF → band structure → DOS/PDOS (plus optional
cutoff and k-point convergence tests), and shows the results interactively: band gap (direct/indirect, VBM/CBM
location), band structure with DOS, projected DOS, SCF convergence, convergence-test recommendations and Löwdin charges.
It can also analyse existing `.out` files.

## Quick start

1. Clone this repository into a Windows folder.
2. Double-click `qe_studio/Start QE Studio.bat`, or in WSL run `bash qe_studio/start.sh`.
   The first start installs Miniforge and QE 7.x from conda-forge (`setup_wsl.sh`).
3. Open http://localhost:8050 and click **Load silicon example** to try it.

See [qe_studio/README.md](qe_studio/README.md) for details.
