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

## Hosting online

The repository includes a `Dockerfile` that runs QE Studio as a website in **public mode**:
- **Blocked on a public server:** reading server folders, opening Explorer and deleting runs.
- **Limited:** each visitor sees only the runs started from their own browser, and jobs have size and time limits.

### Hugging Face Spaces (free: 2 CPUs, 16 GB RAM)

1. Create a free account at https://huggingface.co and a **New Space**: SDK **Docker**, template **Blank**,
   hardware **CPU basic**. Note its id, e.g. `your-name/qe-studio`.
2. Create an access token with **write** permission (Settings → Access Tokens).
3. In this GitHub repository: **Settings → Secrets and variables → Actions**:
   - **Secrets** tab: add `HF_TOKEN` = the token.
   - **Variables** tab: add `HF_SPACE` = the Space id.
4. Run the workflow: **Actions → Deploy to Hugging Face Space → Run workflow**. After that, every push to `main`
   redeploys it. The first build takes about 10 minutes. The app is then at `https://<your-name>-qe-studio.hf.space`.
5. Optional: to restrict access, add a **secret** `QE_STUDIO_PASSWORD` in the Space settings. Visitors then need that
   password (any user name).

Free Spaces sleep after 48 hours without visitors and wake on the next visit. Their disk is temporary, so runs
disappear on restart; download results with **Download all**.

### Any other Docker host (VPS, Render, Railway…)

```bash
docker build -t qe-studio .
docker run -d -p 7860:7860 -e QE_STUDIO_PASSWORD=choose-one -v qe-runs:/home/user/app/runs qe-studio
```

### Settings (environment variables)

| Variable | Public default | Meaning |
|---|---|---|
| `QE_STUDIO_PUBLIC` | `1` in Docker | public mode (restrictions above) |
| `QE_STUDIO_PASSWORD` | – | shared password for the whole site |
| `QE_STUDIO_MAX_ATOMS` | 40 | max atoms per input |
| `QE_STUDIO_MAX_ECUTWFC` | 100 | max ecutwfc (Ry) |
| `QE_STUDIO_MAX_KPOINTS` | 4096 | max points in any k-grid |
| `QE_STUDIO_MAX_STEP_SECONDS` | 900 | time limit per QE program run |
| `QE_STUDIO_MAX_QUEUE` | 5 | max waiting jobs |
| `QE_STUDIO_KEEP_DAYS` | 3 | delete runs older than this |
| `QE_STUDIO_MAX_NP` | all cores | cap on MPI processes |

See [qe_studio/README.md](qe_studio/README.md) for how the app works.
