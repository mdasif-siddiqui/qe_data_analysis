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

**5. NMR (automatic):** NMR needs QE-GIPAW, which is compiled from source. There is nothing to do by hand: the
first time you run a calculation with **NMR** ticked, QE Studio builds it automatically as an extra step
("one-time GIPAW build"), then continues with the NMR calculation.

- **Time:** 30–60 minutes, once per PC; later NMR runs start straight away.
- **Needs:** internet (compilers and source code are downloaded) and about 2.2 GB of disk space. No admin rights.
- **Progress:** the step shows the build messages live; the full logs are in `~/gipaw/*.log`.
- **Skipped when it can't be used:** the build does not start if the input cannot do NMR anyway (smearing instead of
  fixed occupations, or pseudopotentials without GIPAW data).

To build it in advance instead (e.g. before going offline), run once in Ubuntu:

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
| "One-time GIPAW build" step failed | Usually no internet during the build. Check the connection and run the NMR calculation again; the build continues where it stopped. Details: `nmr/gipaw_build.out` in the run folder and `~/gipaw/*.log`. |
| Browser says it can't connect | Wait 2–3 seconds and refresh. If it persists, read the messages in the console window. |
| `Address already in use` | QE Studio is already running in another window; use that one or close it first. |

If a message isn't listed here, copy the full text from the console window.

## Supported materials

QE Studio works with any pw.x `scf.in`, not just the silicon example:

- **Any crystal structure:** every `ibrav`, `celldm` or `A, B, C`, or your own `CELL_PARAMETERS`.
- **Relaxation first** (option **relax first** under SCF: `relax` for atom positions, `vc-relax` for positions and
  cell): a `pw.x` relaxation runs before everything else, with thresholds tight enough for phonons (forces
  < 10⁻⁴ Ry/bohr, pressure < 0.5 kbar). Its final structure is written into the run's `scf.in` (`ibrav = 0`, cell in Å,
  crystal positions; fixed-atom flags kept), and SCF, bands, DOS, phonons, IR & Raman, UV-Vis and NMR all use it.
  The run page reports how the cell, forces and pressure changed. An input that is itself `calculation = 'relax'` or
  `'vc-relax'` is handled the same way. Recommended for phonons and IR & Raman; the app reminds you when they are ticked.
- **Any elements:** pseudopotentials named in `ATOMIC_SPECIES` are downloaded from the official QE library when they exist
  there; otherwise upload your own `.UPF` files together with the input (SSSP, PSlibrary, …).
- **Semiconductors and insulators:** band gap, direct or indirect, with the VBM/CBM positions.
- **Metals** (`occupations = 'smearing'`): detected from bands crossing the Fermi level; plots use E_F = 0.
- **Magnetic materials** (`nspin = 2`): spin-up and spin-down bands and DOS, plus the magnetisation.
- **High-symmetry k-path:** chosen automatically for the lattice type, or entered by hand (e.g. `GXMGR`).
- **Convergence tests** (tick **Convergence tests**): SCF runs at 0.5×, 0.75×, 1×, 1.25×, 1.5× and 2× your ecutwfc
  (same ecutrho/ecutwfc ratio) and on your k-grid scaled by ¼× to 1.5× (your own grid is always one of the points; a
  direction with one k-point, such as a slab's vacuum direction, stays at 1). Energy and pressure differences to the
  largest setting are plotted, and the first value within the threshold (default 1 meV/atom) is recommended. Your own
  lists can be given instead.
- **Phonons** (optional, tick **Phonons**): DFPT with `ph.x` on a q-grid (default 4×4×4), then `q2r.x` and `matdyn.x`.
  - **Dispersion** along the same path as the bands, with the phonon DOS beside it (per element), in cm⁻¹, THz or meV.
  - **Stability check:** warns about imaginary frequencies and large acoustic-sum-rule corrections.
  - **Frequencies at Γ** in a table, before and after the acoustic sum rule.
  - **Harmonic thermodynamics** summed over the phonon modes on a dense q-mesh: zero-point energy, heat capacity
    C_v(T), entropy, free and internal energy from 0 to 1000 K, per mole of atoms.
  - **Phonon DOS** broadening scales with the highest frequency (3 cm⁻¹ for silicon, 9 cm⁻¹ for diamond), so stiff
    materials do not get a spiky DOS.
  - **Cost:** much slower than bands/DOS. Silicon with a 4×4×4 q-grid takes about 20 minutes on 4 cores; for bigger
    cells start with 2×2×2.

- **CIF → scf.in:** drop a `.cif` file (or paste one into the editor).
  - **Structure:** reduced to the primitive cell, with the space group detected (spglib). Positions are symmetrized,
    so rounded CIF values such as `0.33333` become exactly 1/3. A setting suffix on the space-group symbol
    (`F d -3 m :2` for origin choice 2, `R -3 m :R` for rhombohedral axes) is applied when the CIF lists no symmetry
    operations.
  - **Pseudopotentials:** chosen and downloaded from the QE library, either **PBE** (ultrasoft/PAW, PSlibrary) or
    **LDA norm-conserving** (needed for Raman and UV-Vis). Cutoffs come from the pseudopotentials' own recommendations,
    except that norm-conserving B, C, N, O and F get at least 100 Ry (their files suggest too little).
  - **k-grid** from a target spacing (default 0.25 Å⁻¹), with fixed occupations or smearing.
  - **Editing:** the generated input appears in the editor to adjust, and **Save scf.in** downloads it.
  - **Sanity checks:** warns when the built structure does not match the CIF's `_chemical_formula_sum` or its
    `_cell_formula_units_Z` atom count (usually coordinates given for a different origin than the symmetry
    operations), when atoms overlap (< 0.7 Å), and when sites are partially occupied (disorder).
  - **Spin:** *auto* (default) makes the input spin-polarised (`nspin = 2`, ferromagnetic start) when it contains a
    typically magnetic element (Cr, Mn, Fe, Co, Ni, 4f rare earths); *on* / *off* force it. Antiferromagnetic orders
    need a supercell with separate species for up and down sites (not generated automatically).
- **IR & Raman spectra** (tick **IR & Raman**): `ph.x` at Γ with the dielectric tensor and Born effective charges
  (`epsil`), Raman tensors (`lraman`), then `dynmat.x`.
  - **Shown:** ε∞, Born charges, TO and LO frequencies (LO–TO splitting), IR intensities (km/mol) and Raman activities
    (Å⁴/amu) per mode, IR-/Raman-active/silent labels, and broadened spectra. Raman can be shown as intensities for a
    532 nm laser at 300 K.
  - **Cost:** GaAs takes about 5 minutes; wurtzite AlN at 100 Ry about 30 minutes.
- **UV-Vis absorption** (tick **UV-Vis**): an NSCF with many empty bands, then the dielectric function computed
  directly from QE's wavefunctions.
  - **Shown:** ε₁ and ε₂, the absorption coefficient against wavelength (UV and visible ranges), refractive index and
    extinction n and k, reflectivity, the absorption onset, and n at 589 nm.
  - **Scissor shift:** optional, corrects the too-small LDA/GGA gap.
  - **Quality check:** the f-sum rule is reported as a measure of whether enough empty bands were included.
  - **Metals** (smearing, states at E_F): a Drude term for the free carriers is added, with the plasma frequency
    computed from the electron velocities at E_F. The tiles then show the Drude plasma energy, the plasma edge, and
    R, n, k at 589 nm instead of ε∞ and an absorption onset.

- **NMR shielding** (tick **NMR**; the first NMR run builds GIPAW automatically, installation step 5): `pw.x` and `gipaw.x`.
  - **Shown:** isotropic shielding σ, anisotropy, asymmetry η and principal components for every atom, and a simulated
    spectrum per nucleus.
  - **Chemical shifts** use δ = σ_ref − σ. Set σ_ref per element directly, or take it from a reference run plus its
    measured shift (e.g. diamond ¹³C = 35.0 ppm). The app remembers it.
  - **Pseudopotentials:** use **PBE PAW** in the CIF converter. PSlibrary PAW files carry the GIPAW reconstruction data.
  - **Cost:** diamond and SiC take about 1 minute each.

The **Analyse outputs** tab reads existing pw.x, ph.x, dos.x and projwfc.x outputs, matdyn.x `.freq` files
(dispersion) and phonon DOS files (DOS plus thermodynamics) of any material.

## Test: α-quartz from a CIF

[qe_studio/examples/SiO2_alpha-quartz.cif](qe_studio/examples/SiO2_alpha-quartz.cif) is a database-style CIF
(P3₂21, Levien et al. 1980). It was converted with **PBE PAW** and **fixed occupations**, then run with **Band
structure**, **DOS** and **NMR** ticked.

- **Conversion:** 9 atoms (Si₃O₆), space group P3₂21 (No. 154), `O`/`Si.pbe-n-kjpaw_psl.1.0.0.UPF` downloaded
  automatically, cutoffs 50 / 400 Ry, k-grid 6×6×5. The band path is the hexagonal one (Γ–M–K–Γ–A–L–H–A|L–M|K–H).
- **Run time:** about 1.8 hours on 4 cores (SCF 7 min, bands 29 min, NSCF 20 min, NMR 49 min).

| Quantity | QE Studio | Reference |
|---|---|---|
| Band gap | 6.01 eV, indirect (VBM at K, CBM at Γ); direct 6.31 eV | Usual PBE result ≈ 5.7–6 eV; measured ≈ 9 eV (PBE underestimates gaps) |
| ²⁹Si shielding σ | 436.5 ppm (all 3 Si identical) | In line with published GIPAW results for quartz |
| ¹⁷O shielding σ | 217.4 ppm (all 6 O identical); anisotropy 52 ppm, η 0.11 | — |
| GIPAW f-sum rule | −46.2532 | −46.2534 expected ✓ |
| Löwdin charges | Si +1.84, O −0.81 | Partly ionic Si–O bond |
| Warnings | force 0.04 Ry/bohr, pressure 79 kbar | Expected: experimental structure, not relaxed with PBE |

- **²⁹Si chemical shift:** with the 3C-SiC run as reference (σ = 363.88, measured δ = −18.3 ppm), quartz comes out at
  δ ≈ −91 ppm vs −107.4 measured. A reference with very different bonding (Si–C vs Si–O), and computed at other
  cutoffs, is a poor reference. The usual choice is quartz itself: σ_ref = 436.48 − 107.4 = **329.1 ppm** for other
  silicates computed with the same settings.
- **Lesson from the test:** the first version of this CIF had coordinates for another origin than its symmetry
  operations, and ASE silently built a wrong Si₆O₆ cell. The converter now catches this with the composition check.

## Test: ¹⁹F NMR of LiF and NaF from CIFs

[LiF.cif](qe_studio/examples/LiF.cif) and [NaF.cif](qe_studio/examples/NaF.cif) (rock salt, experimental lattices
4.0263 and 4.634 Å) were converted with **PBE PAW** and **fixed occupations** and run with **NMR** ticked. Both used
105 / 840 Ry (the cutoff the converter chose for LiF), so that their ¹⁹F shieldings can be compared.

| | LiF | NaF |
|---|---|---|
| σ(¹⁹F) | 378.00 ppm | 401.68 ppm |
| σ of the cation | ⁷Li 88.27 ppm | ²³Na 533.23 ppm |
| GIPAW f-sum rule | −8.4537 (−8.4533 expected) | −12.7372 (−12.7368 expected) |
| Run time (4 cores) | 3 min | 13 min (CG solver, see below) |

- **Chemical shift:** LiF as the ¹⁹F reference (measured δ = −204.3 ppm vs CFCl₃) gives σ_ref = 173.70 ppm and
  **NaF δ(¹⁹F) = −227.98 ppm, against −224.2 ppm measured**. The 20 ppm difference between the two fluorides is
  reproduced to within 4 ppm. In the app: on the NaF page, pick the LiF run under *Chemical-shift references*, enter
  −204.3 and click **Use**.
- **Fix 1, pseudopotential names:** the converter found no PAW file for Li. PSlibrary names Li `Li.pbe-s-kjpaw`, Mg
  `spnl`, Zn `dnl`; the converter now tries all PSlibrary valence tags (n, dn, spn, spdn, s, sl, nl, dnl, spnl, spdnl,
  sp, spl) and finds Li, Mg and Zn.
- **Fix 2, gipaw.x crash:** for NaF (`Na.pbe-spn-kjpaw`), gipaw.x stopped with "S matrix not positive definite" in
  its Davidson subspace rotation, although the SCF ran fine. QE Studio now detects this and reruns gipaw.x with
  `diagonalization = 'cg'` (slower, but it works); the failed output is kept as `nmr/gipaw_david.out`.

## Test: IR & Raman of wurtzite AlN from a CIF

[qe_studio/examples/AlN_wurtzite.cif](qe_studio/examples/AlN_wurtzite.cif) (P6₃mc, Schulz & Thiemann 1977) was
converted with **LDA norm-conserving** pseudopotentials (`Al`/`N.pz-vbc.UPF`) and **fixed occupations**, then run
with **IR & Raman** ticked. Wurtzite has every kind of mode at Γ: IR + Raman (A1, E1), Raman only (E2), silent (B1),
and a large LO–TO splitting.

- **Settings that matter:** the structure was first relaxed with `vc-relax` in LDA at 100 Ry (a = 3.037 Å, c = 4.848 Å,
  u = 0.383; 2.4% below experiment), and the IR & Raman run used that structure at 100 Ry, k-grid 10×10×6.
  This is now automatic: converting the CIF and ticking **relax first: vc-relax** + **IR & Raman** gives the same cell
  (3.03698 / 4.84823 Å, pressure −136 → −2.4 kbar in 7 SCF cycles, 2 minutes) and the same frequencies to 0.1 cm⁻¹.
- **Run time:** about 40 minutes on 4 cores (ph.x with Raman tensors 27 min).

| Mode | QE Studio (cm⁻¹) | Measured (cm⁻¹) | Label in the app |
|---|---|---|---|
| E2(low) | 251.9 | 248 | Raman |
| B1(low) | 548.9 | — (silent) | silent |
| A1(TO) | 604.4 | 611 | IR + Raman |
| E2(high) | 657.1 | 657 | Raman |
| E1(TO) | 670.1 | 671 | IR + Raman |
| B1(high) | 739.6 | — (silent) | silent |
| E1(LO), q along x | 924.5 | 912 | LO column |

- **Other results:** Born charges Al +2.51 (⊥c) / +2.68 (∥c), N the opposite, summing to 0.001; ε∞ = 4.27 (⊥c) /
  4.51 (∥c), against about 4.7 / 4.8 measured; LDA gap 4.96 eV direct (measured 6.1 eV).
- **What went wrong first:** the converter used the 60 Ry suggested in `N.pz-vbc.UPF`, at the experimental lattice.
  That run gave a pressure of −302 kbar and frequencies about 10% too low (A1(TO) 547, E2(high) 583 cm⁻¹). Between
  60 and 80 Ry the energy still dropped by 82 mRy; 100 Ry is converged (1.5 mRy and 4 kbar more up to 120 Ry). The
  converter now raises norm-conserving first-row elements to 100 Ry and reminds you to relax before IR/Raman.
- **Bug found:** LO frequencies were paired with modes by list position, which put the E1(LO) on a silent B1 mode.
  They are now matched by frequency (GaAs unchanged).

## Test: UV-Vis of diamond from a CIF

[qe_studio/examples/C_diamond.cif](qe_studio/examples/C_diamond.cif) is a minimal CIF without symmetry operations,
with C at (⅛, ⅛, ⅛) and the symbol `F d -3 m :2` (origin choice 2), as many database files are. It was converted with
**LDA norm-conserving** (`C.pz-vbc.UPF`, 100 Ry) and **fixed occupations**, at the experimental lattice
(a = 3.567 Å), then run with **UV-Vis** ticked: 16 bands on a 16×16×16 k-grid, 0.1 eV broadening. Each run takes
about 3 minutes on 4 cores.

| Quantity | No scissor | Scissor +1.7 eV | Measured |
|---|---|---|---|
| Band gap (LDA) | 4.09 eV indirect, 5.62 eV direct | — | 5.47 indirect, 7.3 direct |
| Absorption onset (α = 10⁴ cm⁻¹) | 5.45 eV (227 nm) | 7.14 eV (174 nm) | ≈ 7.3 eV (direct edge) |
| ε∞ | 6.47 | 5.02 | 5.7 |
| n at 589 nm | 2.59 | 2.27 | 2.417 |
| Main ε₂ peak (E2) | 11.8 eV | 13.5 eV | ≈ 12.2 eV |

- **Reading the results:** without a scissor, ε∞ and n are 7–13% too high (the expected overestimate of this
  approximation, see Limitations). The scissor (measured direct gap − LDA direct gap = 1.7 eV) puts the absorption
  edge in the right place but makes ε∞ slightly too low and pushes the main peak too high. The measured values lie
  between the two runs.
- **f-sum rule 77%:** 16 bands are a little few for a spectrum up to 20 eV; more empty bands would raise it.
- **Bug found:** ASE ignored the `:2` suffix, used origin choice 1, and built 16 C atoms instead of 8 (a wrong
  4-atom "diamond", without any warning). The converter now applies the suffix, and a new check compares the atom
  count with `_cell_formula_units_Z`, which also catches wrong settings for one-element crystals.
- **Plot fix:** the absorption plot's wavelength range now follows the onset (diamond: 100–900 nm instead of a
  fixed 150–1500 nm, where its UV absorption was a thin sliver).

## Test: magnetic bcc iron from a CIF

[qe_studio/examples/Fe_bcc.cif](qe_studio/examples/Fe_bcc.cif) (Im-3m, a = 2.8665 Å) was converted with **PBE**
(`Fe.pbe-spn-rrkjus_psl`, 75 / 600 Ry, k-grid 13×13×13), **smearing** and **Spin: auto**, then run with **relax
first: vc-relax**, bands and DOS (about 4.5 minutes on 4 cores).

| | QE Studio | Reference |
|---|---|---|
| Magnetic moment | 2.18 μB/atom (absolute 2.32) | 2.22 μB measured |
| Lattice constant (relaxed, conventional cell) | 2.821 Å | ≈ 2.83 Å typical PBE; 2.8665 Å measured |
| d-band exchange splitting (d-band centres) | 2.21 eV | ≈ 2 eV usually quoted |
| DOS at E_F, spin up / down | 0.84 / 0.21 states/eV | majority much larger than minority |

- **Fix from this test:** the CIF converter never wrote `nspin`, so iron came out as a non-magnetic input. It now has
  the **Spin** option above. Non-magnetic CIFs (AlN, LiF) still convert without spin polarisation.
- **Shown in the app:** the Magnetisation tile, bands split into spin up and spin down, spin-resolved total and
  projected DOS (↓ drawn below the axis), and the relaxation note (pressure −77.9 → −4.8 kbar).
- **Note:** about −5 kbar remains after the relaxation because the final SCF of a vc-relax rebuilds the plane-wave
  basis for the new cell; a higher cutoff reduces this.

## Test: UV-Vis of aluminium (metal) from a CIF

[qe_studio/examples/Al.cif](qe_studio/examples/Al.cif) (fcc, a = 4.0495 Å) was converted with **LDA norm-conserving**
(`Al.pz-vbc`, 50 Ry) and **smearing**, then run with **UV-Vis** ticked: 12 bands on a 16×16×16 k-grid, 0.1 eV
broadening, 15 seconds on 4 cores.

| | QE Studio | Measured / literature |
|---|---|---|
| Plasma edge (bulk plasmon, ε₁ = 0) | 15.64 eV | ≈ 15.0–15.3 eV |
| Drude plasma energy (free carriers) | 12.91 eV | ≈ 12.5–13 eV in other DFT calculations |
| Interband absorption peak | 1.37 eV | ≈ 1.5 eV |
| ε₁ at 2.1 eV (589 nm) | −40.5 | ≈ −50 |
| k / R at 589 nm | 6.37 / 0.99 | ≈ 7 / 0.92 |
| f-sum rule | 98% | 100% ideal |

- **Fix from this test:** UV-Vis had no metal handling: it included only transitions between bands. For aluminium it
  gave ε₁ = +40 to +60 below 1 eV (real Al: large and negative), n(589 nm) = 0.06, an f-sum of 31% and meaningless
  "ε∞" and "onset" tiles. The free-carrier (Drude) term is now computed and added; diamond (an insulator) is unchanged
  (ε∞ = 6.466).
- **Limits:** the Drude damping equals the broadening (0.1 eV here); it is a model parameter, not calculated. With
  it, Al absorbs too little in the visible (R 0.99 instead of 0.92, n too low). The plasma frequency uses a 0.3 eV
  Gaussian at E_F; its convergence with this width and the k-grid was not tested.

## Test: phonons of diamond from a CIF

The same [diamond CIF](qe_studio/examples/C_diamond.cif) (LDA norm-conserving, 100 Ry, k-grid 13×13×13) was first
relaxed with `vc-relax` (a = 3.524 Å, 1.2% below the measured 3.567 Å, typical for LDA; now done by the **relax
first** option), then run with
**Phonons** ticked and the default 4×4×4 q-grid. Diamond is non-polar, so the missing LO–TO splitting in the
dispersion does not matter. Run time: about 35 minutes on 4 cores (ph.x 31 min for 8 q-points).

| Frequencies (cm⁻¹) | QE Studio | Measured |
|---|---|---|
| Γ optical (Raman line) | 1335.3 | 1332 |
| X: TA / TO / LA = LO | 795.5 / 1105.9 / 1229.3 | 807 / 1072 / 1185 |
| L: TA / LA / TO / LO | 555.8 / 1083.8 / 1241.4 / 1275.4 | 552 / 1035 / 1210 / 1242 |
| Highest frequency (overbending, away from Γ) | 1355 | seen in experiment |
| Phonon DOS main peak | 1250 | ≈ 1250 |

| Thermodynamics | QE Studio | Measured / literature |
|---|---|---|
| C_v at 300 K | 6.02 J/(mol·K) | ≈ 6.1 |
| S at 300 K | 2.32 J/(mol·K) | 2.38 |
| Zero-point energy | 185 meV/atom | ≈ 180 |
| C_v at 1000 K | 20.9 J/(mol·K) | ≈ 21.3 |

- **Accuracy:** Γ within 0.2% and acoustic modes within 1.5%; optical modes at X and L 3–5% high, which is in the
  usual range for LDA with the old `C.pz-vbc` pseudopotential. Stable (no imaginary frequencies); the acoustic modes
  at Γ were 16 cm⁻¹ before the acoustic sum rule, which the app applies.
- **Fix from this test:** the phonon DOS used a fixed 3 cm⁻¹ broadening, made for silicon (up to 520 cm⁻¹). On the
  same q-mesh diamond's DOS (up to 1355 cm⁻¹) was very spiky. The broadening is now (highest frequency)/150, at
  least 3 cm⁻¹: 15 times smoother for diamond at the same cost, silicon essentially unchanged.

## Test: convergence tests on AlN from a CIF

The [AlN CIF](qe_studio/examples/AlN_wurtzite.cif) was converted with **LDA norm-conserving** (100 / 400 Ry, k-grid
10×10×6) and run with **Convergence tests** ticked and the default settings (about 8 minutes on 4 cores). The
energies agree to all digits with the SCF runs made by hand for the IR & Raman test (e.g. −48.71566948 Ry at 100 Ry).

| ecutwfc (Ry) | 50 | 75 | 100 | 125 | 150 | 200 |
|---|---|---|---|---|---|---|
| ΔE to 200 Ry (meV/atom) | 768 | 83 | 12.1 | 6.4 | 2.9 | 0 |
| ΔP to 200 Ry (kbar) | −303 | −58 | −5.1 | −1.4 | −1.5 | 0 |

| k-grid | 2×2×2 | 5×5×3 | 8×8×4 | 10×10×6 | 12×12×8 | 15×15×9 |
|---|---|---|---|---|---|---|
| ΔE to 15×15×9 (meV/atom) | 305 | 1.95 | 0.014 | 0.005 | −0.003 | 0 |

- **Cutoff:** the app reports that the total energy is not converged to 1 meV/atom up to 200 Ry, which is correct for
  the old `N.pz-vbc` pseudopotential. Pressure (and with it structures and frequencies) is converged much earlier, to
  about 1.5 kbar from 125 Ry; only absolute total energies need the very high cutoff.
- **k-points:** recommended 8×8×4 (insulator, converges fast).
- **Fix from this test:** the k test used to rebuild its grids from a k-point density, so its "1×" point was 10×10×5
  instead of the input's 10×10×6. It now scales the input's own grid, as in the table above.

## Limitations

- **Testing so far:** silicon (semiconductor) and aluminium (metal) have been run end to end; both are built-in
  examples. α-quartz and LiF/NaF (NMR), wurtzite AlN (IR & Raman, convergence tests) and diamond (UV-Vis, phonons), aluminium (UV-Vis of a metal) and bcc iron (ferromagnet) were run end to end from CIFs
  (see above). Input reading was also checked for a hexagonal 2D slab and for GaAs with `ibrav = 0`. Antiferromagnets
  and other materials should work but are not yet tested; check the first results critically.
- **Metals need dense k-grids.** The aluminium example's 12×12×12 grid is deliberately not converged: the app warns
  that E_F moves by +130 meV on the denser NSCF grid and that the convergence tests did not converge.
- **2D materials / slabs:** the automatic k-path includes points along the vacuum direction. Enter the in-plane
  path yourself, e.g. `GMKG` for a hexagonal layer.
- **Relaxation:** the relaxed structure is used as it is; there is no symmetry re-detection or re-standardization
  afterwards. A relaxation that does not converge within 100 steps is reported, and its last structure is used.
- **Spin-orbit / non-collinear** (`noncolin`, `lspinorb`): basic support only, the least tested part.
- **Large systems:** on a 4-core laptop, more than about 20–30 atoms or high cutoffs take from many minutes to hours.
- **Phonons:**
  - **Physics not included:** the phonon dispersion has no LO–TO splitting, so optical branches near Γ of polar
    materials (e.g. GaAs, NaCl) are not fully correct; the **IR & Raman** option does give the LO and TO frequencies
    at Γ. No anharmonic effects or thermal expansion.
  - **Inputs:** phonons need a relaxed structure (zero forces) and tight convergence; the app runs its own SCF with
    `conv_thr` ≤ 1e-12.
  - **Relax first:** choose **relax first: vc-relax** (or give an input that is already relaxed with the same
    pseudopotentials and cutoff).
  - **Testing:** silicon and diamond (from a CIF, see above).
- **IR & Raman:**
  - **Insulators only** (`occupations = 'fixed'`).
  - **Raman** also needs LDA with norm-conserving pseudopotentials (a ph.x restriction); otherwise only IR is calculated.
  - **Spectra:** harmonic, at Γ only, with no two-phonon or anharmonic features.
  - **Testing:** done with GaAs (TO 258 / LO 274 cm⁻¹ in LDA vs 268 / 292 measured) and wurtzite AlN from a CIF
    (within 1–2% of measured, see above).
  - **Relax first:** choose **relax first: vc-relax**; at the experimental lattice frequencies can be off by 10%.
- **UV-Vis:**
  - **Approximation:** independent-particle, so no excitons and no local-field effects.
  - **Accuracy:** the nonlocal-pseudopotential commutator is left out, which overestimates ε by about 10–20% (silicon:
    ε∞ = 16.2 vs about 14 in full LDA calculations, 11.7 measured; diamond: 6.47 vs 5.7 measured).
  - **Band gaps:** LDA/GGA red-shift the absorption edge; use the scissor shift.
  - **Requirements:** norm-conserving pseudopotentials and QE's HDF5 wavefunctions (the conda-forge build has them).
  - **Why not epsilon.x:** QE's own `epsilon.x` is not used; in this build it gives k-grid-dependent spectra that are
    several times too large.
  - **Metals:** the Drude damping is not calculated (it equals the broadening), so visible absorption and n of metals
    are only qualitative; the plasma edge and ε₁ are more reliable. Tested with aluminium.
- **NMR:**
  - **Insulators only** (`occupations = 'fixed'`). Knight shifts of metals are not calculated.
  - **Pseudopotentials:** they must carry GIPAW data. PSlibrary **PAW** files work. PSlibrary ultrasoft files contain
    GIPAW sections with zero radii that `gipaw.x` rejects, and norm-conserving `pz-vbc`/`bhs` files have none; the app
    checks this and explains.
  - **Accuracy:** chemical shifts are only as good as the reference. Use a reference computed with the same
    pseudopotentials, cutoff and k-grid. Tested: 3C-SiC ¹³C δ = 20.2 ppm with diamond as reference, vs 18.4 measured;
    α-quartz ²⁹Si σ = 436.5 ppm (but δ ≈ −91 vs −107.4 ppm when SiC is used as the reference); NaF ¹⁹F δ = −228.0
    vs −224.2 ppm with LiF as the reference.
  - **Solver fallback:** if gipaw.x stops with "S matrix not positive definite", it is rerun automatically with the
    CG solver, which takes several times longer (NaF: 13 instead of about 3 minutes).
  - **Build patches:** `build_gipaw.sh` fixes three malformed Fortran format strings in QE-GIPAW 7.6 that make recent
    gfortran stop at the end of a run.
- **Other QE programs** (`cp.x`, `neb.x`, `pp.x`, …) are not supported.
- **`ATOMIC_POSITIONS crystal_sg`** is not supported by the structure viewer.
- **Spin-polarised or non-collinear** UV-Vis, and NMR of magnetic systems, are untested.

See [qe_studio/README.md](qe_studio/README.md) for how the app works.
