#!/usr/bin/env bash
# Build QE-GIPAW (NMR shielding) inside WSL, without admin rights.
# Creates a conda env "qebuild" with compilers + MPI + FFTW + OpenBLAS, compiles QE's pw.x and gipaw.x,
# and installs both into ~/gipaw/bin. Takes about 30-60 minutes. Safe to re-run: finished steps are skipped.
set -euo pipefail

QE_VER="7.6"                     # must match a qe-gipaw tag
MF="$HOME/miniforge3"
ROOT="$HOME/gipaw"
SRC="$ROOT/src"
mkdir -p "$SRC" "$ROOT/bin"

if ! "$MF/bin/conda" env list | grep -q '^qebuild '; then
  echo ">>> Creating build env 'qebuild' (compilers, MPI, FFTW, OpenBLAS)"
  "$MF/bin/conda" create -n qebuild -c conda-forge -y compilers make openmpi fftw openblas git
fi
ENV="$MF/envs/qebuild"
export PATH="$ENV/bin:$PATH"
export LD_LIBRARY_PATH="$ENV/lib:${LD_LIBRARY_PATH:-}"

cd "$SRC"
if [ ! -d "q-e-qe-$QE_VER" ]; then
  echo ">>> Downloading Quantum ESPRESSO $QE_VER source"
  curl -fsSL "https://gitlab.com/QEF/q-e/-/archive/qe-$QE_VER/q-e-qe-$QE_VER.tar.gz" | tar xz
fi
QE="$SRC/q-e-qe-$QE_VER"
cd "$QE"
if [ ! -x bin/pw.x ]; then
  echo ">>> Configuring QE"
  make veryclean > /dev/null 2>&1 || true
  # FFTW is auto-detected (giving FFT_LIBS by hand skips the -D__FFTW3 flag and the build fails)
  ./configure MPIF90=mpif90 CC=mpicc F90=gfortran \
      BLAS_LIBS="-L$ENV/lib -lopenblas" LAPACK_LIBS="-L$ENV/lib -lopenblas" \
      CPPFLAGS="-I$ENV/include" LIBDIRS="$ENV/lib" \
      LDFLAGS="-L$ENV/lib -Wl,-rpath,$ENV/lib" --with-scalapack=no > "$ROOT/configure_qe.log" 2>&1
  if ! grep -q "__FFTW3" make.inc; then
    echo "ERROR: configure did not find FFTW3 (see $ROOT/configure_qe.log)"; exit 1
  fi
  echo ">>> Compiling pw.x (about 15-30 minutes)"
  make -j"$(nproc --all)" pw > "$ROOT/make_pw.log" 2>&1
fi
test -x bin/pw.x

if [ ! -d "$QE/qe-gipaw" ]; then
  echo ">>> Downloading QE-GIPAW $QE_VER"
  git clone -q --depth 1 --branch "$QE_VER" https://github.com/dceresoli/qe-gipaw.git "$QE/qe-gipaw"
fi
cd "$QE/qe-gipaw"
# Some GIPAW format strings lack a comma before a closing string, e.g. i3")" instead of i3,")".
# Recent gfortran stops at run time with "Missing comma between descriptors"; add the commas.
if grep -qE '(i[0-9]+)"\)"' src/gipaw_routines.f90; then
  sed -i -E 's/(i[0-9]+)"\)"/\1,")"/g' src/gipaw_routines.f90
  rm -f src/gipaw.x src/gipaw_routines.o
fi
for f in suscept_crystal core_relax; do                  # same problem with '' quoting: F10.4''): and I1''S
  if grep -qE "([FI][0-9]+(\.[0-9]+)?)''" "src/$f.f90"; then
    sed -i -E "s/([FI][0-9]+(\.[0-9]+)?)''/\1,''/g" "src/$f.f90"
    rm -f src/gipaw.x "src/$f.o"
  fi
done
if [ ! -x src/gipaw.x ]; then
  echo ">>> Compiling gipaw.x"
  # GIPAW links QE's linear-response library; build it first (a parallel make otherwise races for it)
  (cd "$QE" && make -j"$(nproc --all)" lrmods > "$ROOT/make_lrmods.log" 2>&1)
  ./configure --with-qe-source="$QE" > "$ROOT/configure_gipaw.log" 2>&1
  make -j"$(nproc --all)" > "$ROOT/make_gipaw.log" 2>&1
fi
test -x src/gipaw.x

cp -f "$QE/bin/pw.x" "$ROOT/bin/pw.x"
cp -f "$QE/qe-gipaw/src/gipaw.x" "$ROOT/bin/gipaw.x"
cat > "$ROOT/env.sh" <<EOF
export PATH="$ROOT/bin:$ENV/bin:\$PATH"
export LD_LIBRARY_PATH="$ENV/lib:\${LD_LIBRARY_PATH:-}"
EOF
echo ">>> Done: $ROOT/bin/pw.x and $ROOT/bin/gipaw.x"
echo "" | OMP_NUM_THREADS=1 "$ROOT/bin/pw.x" 2>/dev/null | grep -m1 "Program PWSCF" || true
