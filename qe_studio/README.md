# QE Studio

A local web app for Quantum ESPRESSO: give it a pw.x `scf.in`, and it checks the input, runs the
calculations in WSL and shows the results as interactive plots.

## Start

- **Windows:** double-click `Start QE Studio.bat` (opens http://localhost:8050).
- **WSL terminal:** `bash start.sh`, then open http://localhost:8050.

The first start runs `setup_wsl.sh`, which installs Miniforge and a conda env `qe` with QE 7.x from
conda-forge. It never uses apt's `quantum-espresso`, which is broken on Ubuntu 24.04.

## What it does

**New calculation**: drop an `scf.in` (and optionally `.UPF` files) or paste it into the editor.
- Checks the input: namelists, cards, ecutrho ≥ 4·ecutwfc, missing cards, post-processing input given by mistake.
  It also fixes CRLF line endings and a missing final newline.
- Shows the structure (any `ibrav`, `celldm` or `A,B,C`, `CELL_PARAMETERS`) in 3D, with lattice type, k-point spacing and the suggested high-symmetry path.
- Finds pseudopotentials: uploaded files, then `pseudo_library/`, then `$ESPRESSO_PSEUDO`. Missing ones are
  downloaded from pseudopotentials.quantum-espresso.org, checked to be real UPF files (not an HTML error page) and cached.
- Runs, in order: SCF → bands (pw.x along the k-path, crystal_b) → bands.x → NSCF on a dense grid → dos.x → projwfc.x,
  plus optional ecutwfc and k-point convergence tests.

**Results page** for each run:
- Band gap (direct/indirect) with VBM/CBM positions (for example "85% of the way Γ→X"), total energy, HOMO/LUMO or E_F, pressure, forces, magnetisation.
- Band structure with DOS/PDOS on a shared energy axis (zero at the VBM for insulators, at E_F for metals).
- Projected DOS by element and orbital, SCF convergence, convergence-test plots with the recommended value, Löwdin charges.
- Checks: every output ends with JOB DONE, convergence NOT achieved, unrelaxed forces, non-zero pressure.
- All files are downloadable, or open the run folder in Explorer. Plots export to PNG with the camera button.

**Analyse outputs**: drop existing `.out` files (and the bands `.in` for k-labels), or give a folder path
(`C:\...` or `/mnt/c/...`). It produces a summary table with CSV export, plots of SCF convergence, an energy comparison, relaxation
trajectories, ecut/k convergence across files and band structures.

## Layout

```
app.py              Flask server (127.0.0.1:8050)
qestudio/qein.py    input parsing/editing, ibrav → cell (port of QE latgen), k-paths (ASE)
qestudio/qeout.py   pw.x / dos.x / projwfc.x parsers, band-gap analysis
qestudio/jobs.py    job queue + pipeline, MPI launch fallback, convergence tests, analysis
qestudio/pseudos.py pseudopotential lookup / download / validation
static/             UI (HTML/CSS/JS + Plotly)
runs/<id>/          one folder per run: inputs, outputs, results.json, summary.csv
```

## Notes

- MPI: `mpirun -np <physical cores>` with `OMP_NUM_THREADS=1`. If launching fails, it retries with CMA disabled, then serially.
  Using hyper-threads (8 ranks on this 4-core i7) was ~10× slower in testing, so the default is the physical core count.
- QE's `outdir` (wavefunctions) goes to `~/.cache/qe_studio/<run>` on the Linux disk, which is much faster than `/mnt/c`, and is deleted when
  the run ends. Inputs and outputs stay in `runs/`.
- Jobs run one at a time. Cancel stops the running program.
- relax/vc-relax inputs are accepted, but the bands/DOS inputs are generated from your input file's structure.
  For band structures of a relaxed geometry, it is safest to put the final coordinates into a new scf input and run that.
