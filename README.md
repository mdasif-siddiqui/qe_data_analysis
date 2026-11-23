# QE Data Analysis

**QE Studio** is a local web app for running and analysing [Quantum ESPRESSO](https://www.quantum-espresso.org/)
calculations on Windows via WSL Ubuntu.

Give it a pw.x `scf.in` and it validates the input, runs SCF → band structure → DOS/PDOS (plus optional
cutoff and k-point convergence tests), and shows the results interactively: band gap (direct/indirect, VBM/CBM
location), band structure with DOS, projected DOS, SCF convergence, convergence-test recommendations and Löwdin charges.
It can also analyse existing `.out` files.

## Installation (Windows 10/11)

**1. Install WSL with Ubuntu** (once, if the PC doesn't have it). In PowerShell as Administrator:

```powershell
wsl --install -d Ubuntu
```

Restart when asked, then open Ubuntu once and create a Linux username and password.

**2. Get the project.** With Git for Windows installed:

```powershell
git clone https://github.com/mdasif-siddiqui/qe_data_analysis.git C:\Users\$env:USERNAME\qe_data_analysis
```

Without Git, use **Code → Download ZIP** on GitHub and unzip it.

**3. Start it.** Double-click `qe_data_analysis\qe_studio\Start QE Studio.bat`.

- **First start:** `setup_wsl.sh` installs Miniforge plus QE 7.x and the Python packages from conda-forge inside WSL.
  This needs internet, about 2 GB of disk, and takes about 5–10 minutes.
- **Every start:** the browser opens http://localhost:8050. If the page loads before the server is ready, refresh it.
  Keep the console window open while you use the app; closing it stops the server.

**4. Test it.** Click **Load silicon example**, then **Run calculation**. You should get a band gap of about 0.50 eV (indirect).

### What happens automatically

- **MPI processes:** set to the PC's physical core count.
- **Pseudopotentials:** downloaded from the QE library on first use and cached in `qe_studio/pseudo_library/`.
- **Runs:** saved in `qe_studio/runs/`. Neither runs nor pseudopotentials are in the repository; to move old results
  to another PC, copy the `qe_studio\runs\` folder.

### Other systems

- **Native Linux:** clone the repo and run `bash qe_studio/start.sh`. Everything works except the "Open folder"
  button, which only opens Windows Explorer.
- **macOS:** not supported as-is. `setup_wsl.sh` installs the Linux build of Miniforge.

### Keeping several PCs in sync

After changing code on one PC, commit and push there, then run `git pull` in `qe_data_analysis` on the other.

## Supported materials

QE Studio works with any pw.x `scf.in`, not just the silicon example:

- **Any crystal structure:** every `ibrav`, `celldm` or `A, B, C`, or your own `CELL_PARAMETERS`.
- **Any elements:** pseudopotentials named in `ATOMIC_SPECIES` are downloaded from the official QE library when they exist
  there; otherwise upload your own `.UPF` files together with the input (SSSP, PSlibrary, …).
- **Semiconductors and insulators:** band gap, direct or indirect, with the VBM/CBM positions.
- **Metals** (`occupations = 'smearing'`): detected from bands crossing the Fermi level; plots use E_F = 0.
- **Magnetic materials** (`nspin = 2`): spin-up and spin-down bands and DOS, plus the magnetisation.
- **High-symmetry k-path:** chosen automatically for the lattice type, or entered by hand (e.g. `GXMGR`).

The **Analyse outputs** tab reads existing pw.x, dos.x and projwfc.x outputs of any material.

## Limitations

- **Testing so far:** only silicon has been run end to end. Input reading was also checked for a hexagonal 2D slab
  and for GaAs with `ibrav = 0`. Metals, magnetic systems and other materials should work but are not yet tested;
  check the first results critically.
- **2D materials / slabs:** the automatic k-path includes points along the vacuum direction. Enter the in-plane
  path yourself, e.g. `GMKG` for a hexagonal layer.
- **relax / vc-relax inputs:** they run, but the band structure and DOS inputs use the structure from your input file,
  not the relaxed one. For bands of the relaxed structure, put the final coordinates into a new scf input.
- **Spin-orbit / non-collinear** (`noncolin`, `lspinorb`): basic support only, the least tested part.
- **Large systems:** on a 4-core laptop, more than about 20–30 atoms or high cutoffs take from many minutes to hours.
- **Other QE programs** (`ph.x`, `cp.x`, `neb.x`, …) are not supported.
- **`ATOMIC_POSITIONS crystal_sg`** is not supported by the structure viewer.

See [qe_studio/README.md](qe_studio/README.md) for how the app works.
