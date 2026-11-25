# QE Data Analysis

**QE Studio** is a local web app for running and analysing [Quantum ESPRESSO](https://www.quantum-espresso.org/)
calculations on Windows via WSL Ubuntu.

Give it a pw.x `scf.in` (or a `.cif` file, which it converts) and it validates the input, runs SCF → band structure →
DOS/PDOS, plus optional phonons, IR and Raman spectra, UV-Vis absorption, NMR shielding and cutoff/k-point convergence
tests. It shows
the results interactively:
- **Electronic:** band gap (direct/indirect, VBM/CBM location), band structure with DOS, projected DOS, Löwdin charges.
- **Vibrational:** phonon dispersion and DOS, heat capacity and other harmonic thermodynamics, IR and Raman spectra.
- **Optical:** dielectric function, absorption coefficient, refractive index.
- **NMR:** shielding tensors and chemical shifts for every atom.
- **Quality:** SCF convergence and convergence-test recommendations.

It can also analyse existing output files.

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

**5. Optional, only for NMR:** in Ubuntu, run once:

```bash
bash /mnt/c/Users/<your-name>/qe_data_analysis/qe_studio/build_gipaw.sh
```

(The project folder on the Windows drive is under `/mnt/c/...` inside Ubuntu.)

This compiles QE-GIPAW (`gipaw.x`, not available in conda-forge) into `~/gipaw`, in its own conda environment.
It needs no admin rights, takes about 30–60 minutes and about 2 GB of disk. Restart QE Studio afterwards.

## Daily use

After the first start, nothing is installed again:

1. **Double-click `qe_studio\Start QE Studio.bat`.** The server starts in a few seconds and the browser opens
   http://localhost:8050 by itself. If the page says it can't connect, wait 2–3 seconds and refresh.
2. **Use your own input.** Drop your `scf.in` or a `.cif` file (and any `.UPF` files) into the page, choose what to calculate and click
   **Run calculation**. The built-in examples are only for checking the installation.
3. **Look at old results anytime.** The **Runs** tab lists every calculation, saved in `qe_studio\runs\`.
4. **When you're done, close the console window.** That stops the app. Don't close it while a calculation is
   running, or the calculation stops too.

In short: **double-click → work in the browser → close the window.**

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

## Troubleshooting

Start with `git pull` (or download the ZIP again): the start scripts fix several of these problems by themselves.

| Message | Cause and fix |
|---|---|
| `python: not found` or `No module named flask` | QE or its Python packages are missing, for example because a `qe` conda env made by hand (`conda create -n qe qe`) has no Python. The current `start.sh` detects this and repairs the env automatically (needs internet). |
| `ERROR: could not start Ubuntu in WSL` | Ubuntu isn't installed or isn't the default. In PowerShell, run `wsl -l -v`: Ubuntu must have the `*`. Fix with `wsl --set-default Ubuntu` (use the exact name shown, e.g. `Ubuntu-22.04`). |
| `$'\r': command not found` | The `.sh` files got Windows line endings (copied or edited in Windows). `Start QE Studio.bat` now fixes this on every start. By hand: `sed -i 's/\r$//' qe_studio/*.sh` in Ubuntu. |
| `ERROR: curl is missing` | In Ubuntu, run `sudo apt update && sudo apt install -y curl`, then start again. |
| `Setup failed` | Usually no internet during the first start. Check the connection and start again; setup continues where it stopped. |
| NMR option greyed out | GIPAW isn't built on this PC: run `build_gipaw.sh` once in Ubuntu (installation step 5), then restart QE Studio. If it fails, read `~/gipaw/*.log`. |
| Browser says it can't connect | Wait 2–3 seconds and refresh. If it persists, read the messages in the console window. |
| `Address already in use` | QE Studio is already running in another window; use that one or close it first. |

If a message isn't listed here, copy the full text from the console window.

## Supported materials

QE Studio works with any pw.x `scf.in`, not just the silicon example:

- **Any crystal structure:** every `ibrav`, `celldm` or `A, B, C`, or your own `CELL_PARAMETERS`.
- **Any elements:** pseudopotentials named in `ATOMIC_SPECIES` are downloaded from the official QE library when they exist
  there; otherwise upload your own `.UPF` files together with the input (SSSP, PSlibrary, …).
- **Semiconductors and insulators:** band gap, direct or indirect, with the VBM/CBM positions.
- **Metals** (`occupations = 'smearing'`): detected from bands crossing the Fermi level; plots use E_F = 0.
- **Magnetic materials** (`nspin = 2`): spin-up and spin-down bands and DOS, plus the magnetisation.
- **High-symmetry k-path:** chosen automatically for the lattice type, or entered by hand (e.g. `GXMGR`).
- **Phonons** (optional, tick **Phonons**): DFPT with `ph.x` on a q-grid (default 4×4×4), then `q2r.x` and `matdyn.x`.
  - **Dispersion** along the same path as the bands, with the phonon DOS beside it (per element), in cm⁻¹, THz or meV.
  - **Stability check:** warns about imaginary frequencies and large acoustic-sum-rule corrections.
  - **Frequencies at Γ** in a table, before and after the acoustic sum rule.
  - **Harmonic thermodynamics** from the phonon DOS: zero-point energy, heat capacity C_v(T), entropy, free and internal
    energy from 0 to 1000 K, per mole of atoms.
  - **Cost:** much slower than bands/DOS. Silicon with a 4×4×4 q-grid takes about 20 minutes on 4 cores; for bigger
    cells start with 2×2×2.

- **CIF → scf.in:** drop a `.cif` file (or paste one into the editor).
  - **Structure:** reduced to the primitive cell, with the space group detected (spglib).
  - **Pseudopotentials:** chosen and downloaded from the QE library, either **PBE** (ultrasoft/PAW, PSlibrary) or
    **LDA norm-conserving** (needed for Raman and UV-Vis). Cutoffs come from the pseudopotentials' own recommendations.
  - **k-grid** from a target spacing (default 0.25 Å⁻¹), with fixed occupations or smearing.
  - **Editing:** the generated input appears in the editor to adjust, and **Save scf.in** downloads it.
- **IR & Raman spectra** (tick **IR & Raman**): `ph.x` at Γ with the dielectric tensor and Born effective charges
  (`epsil`), Raman tensors (`lraman`), then `dynmat.x`.
  - **Shown:** ε∞, Born charges, TO and LO frequencies (LO–TO splitting), IR intensities (km/mol) and Raman activities
    (Å⁴/amu) per mode, IR-/Raman-active/silent labels, and broadened spectra. Raman can be shown as intensities for a
    532 nm laser at 300 K.
  - **Cost:** GaAs takes about 5 minutes.
- **UV-Vis absorption** (tick **UV-Vis**): an NSCF with many empty bands, then the dielectric function computed
  directly from QE's wavefunctions.
  - **Shown:** ε₁ and ε₂, the absorption coefficient against wavelength (UV and visible ranges), refractive index and
    extinction n and k, reflectivity, the absorption onset, and n at 589 nm.
  - **Scissor shift:** optional, corrects the too-small LDA/GGA gap.
  - **Quality check:** the f-sum rule is reported as a measure of whether enough empty bands were included.

- **NMR shielding** (tick **NMR**; needs the one-time GIPAW build, installation step 5): `pw.x` and `gipaw.x`.
  - **Shown:** isotropic shielding σ, anisotropy, asymmetry η and principal components for every atom, and a simulated
    spectrum per nucleus.
  - **Chemical shifts** use δ = σ_ref − σ. Set σ_ref per element directly, or take it from a reference run plus its
    measured shift (e.g. diamond ¹³C = 35.0 ppm). The app remembers it.
  - **Pseudopotentials:** use **PBE PAW** in the CIF converter. PSlibrary PAW files carry the GIPAW reconstruction data.
  - **Cost:** diamond and SiC take about 1 minute each.

The **Analyse outputs** tab reads existing pw.x, ph.x, dos.x and projwfc.x outputs, matdyn.x `.freq` files
(dispersion) and phonon DOS files (DOS plus thermodynamics) of any material.

## Limitations

- **Testing so far:** silicon (semiconductor) and aluminium (metal) have been run end to end; both are built-in
  examples. Input reading was also checked for a hexagonal 2D slab and for GaAs with `ibrav = 0`. Magnetic systems
  and other materials should work but are not yet tested; check the first results critically.
- **Metals need dense k-grids.** The aluminium example's 12×12×12 grid is deliberately not converged: the app warns
  that E_F moves by +130 meV on the denser NSCF grid and that the convergence tests did not converge.
- **2D materials / slabs:** the automatic k-path includes points along the vacuum direction. Enter the in-plane
  path yourself, e.g. `GMKG` for a hexagonal layer.
- **relax / vc-relax inputs:** they run, but the band structure and DOS inputs use the structure from your input file,
  not the relaxed one. For bands of the relaxed structure, put the final coordinates into a new scf input.
- **Spin-orbit / non-collinear** (`noncolin`, `lspinorb`): basic support only, the least tested part.
- **Large systems:** on a 4-core laptop, more than about 20–30 atoms or high cutoffs take from many minutes to hours.
- **Phonons:**
  - **Physics not included:** the phonon dispersion has no LO–TO splitting, so optical branches near Γ of polar
    materials (e.g. GaAs, NaCl) are not fully correct; the **IR & Raman** option does give the LO and TO frequencies
    at Γ. No anharmonic effects or thermal expansion.
  - **Inputs:** phonons need a relaxed structure (zero forces) and tight convergence; the app runs its own SCF with
    `conv_thr` ≤ 1e-12.
  - **Testing:** only silicon has been tested with phonons so far.
- **IR & Raman:**
  - **Insulators only** (`occupations = 'fixed'`).
  - **Raman** also needs LDA with norm-conserving pseudopotentials (a ph.x restriction); otherwise only IR is calculated.
  - **Spectra:** harmonic, at Γ only, with no two-phonon or anharmonic features.
  - **Testing:** done with GaAs (TO 258 / LO 274 cm⁻¹ in LDA vs 268 / 292 measured).
- **UV-Vis:**
  - **Approximation:** independent-particle, so no excitons and no local-field effects.
  - **Accuracy:** the nonlocal-pseudopotential commutator is left out, which overestimates ε by about 10–20% (silicon:
    ε∞ = 16.2 vs about 14 in full LDA calculations, 11.7 measured).
  - **Band gaps:** LDA/GGA red-shift the absorption edge; use the scissor shift.
  - **Requirements:** norm-conserving pseudopotentials and QE's HDF5 wavefunctions (the conda-forge build has them).
  - **Why not epsilon.x:** QE's own `epsilon.x` is not used; in this build it gives k-grid-dependent spectra that are
    several times too large.
- **NMR:**
  - **Insulators only** (`occupations = 'fixed'`). Knight shifts of metals are not calculated.
  - **Pseudopotentials:** they must carry GIPAW data. PSlibrary **PAW** files work. PSlibrary ultrasoft files contain
    GIPAW sections with zero radii that `gipaw.x` rejects, and norm-conserving `pz-vbc`/`bhs` files have none; the app
    checks this and explains.
  - **Accuracy:** chemical shifts are only as good as the reference. Use a reference computed with the same
    pseudopotentials, cutoff and k-grid. Tested: 3C-SiC ¹³C δ = 20.2 ppm with diamond as reference, vs 18.4 measured.
  - **Build patches:** `build_gipaw.sh` fixes three malformed Fortran format strings in QE-GIPAW 7.6 that make recent
    gfortran stop at the end of a run.
- **Other QE programs** (`cp.x`, `neb.x`, `pp.x`, …) are not supported.
- **`ATOMIC_POSITIONS crystal_sg`** is not supported by the structure viewer.
- **Spin-polarised or non-collinear** UV-Vis, and NMR of magnetic systems, are untested.

See [qe_studio/README.md](qe_studio/README.md) for how the app works.
