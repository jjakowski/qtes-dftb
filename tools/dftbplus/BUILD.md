# Building the DFTB+ library for qtes-dftb

DFTB+ is maintained externally and built **outside this repository**. These
are the exact commands used for milestone 1 (2026-11-17). See
`docs/DFTBPLUS_AUDIT.md` for the background.

| Item | Value |
|---|---|
| DFTB+ version | 25.1 (tag `25.1`, commit `8802aff6`), library API 0.4.0 |
| Source | `/home/j2c/Development/DFTB/dftbplus-25.1/` (exported from the checkout `/home/j2c/Development/DFTB/dftbplus`, which stays untouched) |
| Build directory | `/home/j2c/Development/DFTB/dftbplus-25.1/build-qtes/` |
| Install prefix | `/home/j2c/Development/DFTB/dftbplus-25.1/install-qtes/` |
| Toolchain | conda env `rmg-dev`: gfortran/gcc 15.2.0, CMake 4.1.2, OpenBLAS |
| Options | `WITH_API=ON`, `WITH_MPI=OFF`, `WITH_OMP=OFF`, static library, Release |

## Commands

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate rmg-dev

# 1. Export the release (reads the checkout, changes nothing in it)
B=/home/j2c/Development/DFTB/dftbplus-25.1
mkdir -p $B
cd /home/j2c/Development/DFTB/dftbplus
git archive 25.1 | tar -x -C $B

# 2. Configure (about 1 s; no downloads needed for this configuration)
cd $B
FC=gfortran CC=gcc CMAKE_PREFIX_PATH=$CONDA_PREFIX cmake -S . -B build-qtes \
  -DWITH_API=ON -DWITH_MPI=OFF -DWITH_OMP=OFF -DBUILD_SHARED_LIBS=OFF \
  -DBLAS_LIBRARY=openblas -DLAPACK_LIBRARY=openblas \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$B/install-qtes > cmake-qtes.log 2>&1

# 3. Build (about 1 min with -j 32) and install
cmake --build build-qtes -j 32 > build-qtes.log 2>&1
cmake --install build-qtes > install-qtes.log 2>&1
```

The logs `cmake-qtes.log`, `build-qtes.log` and `install-qtes.log` are in
`$B`.

## What gets installed

| Path (under `install-qtes/`) | Content |
|---|---|
| `lib/libdftbplus.a` | the library (static) |
| `include/dftbplus/modfiles/` | 299 Fortran module files, incl. `dftbplus.mod` (the API) |
| `include/dftbplus.h` | C API header |
| `lib/pkgconfig/dftbplus.pc` | compile and link flags |
| `lib/cmake/dftbplus/` | CMake package config |
| `bin/dftb+` | the standalone program (same build), used for comparisons |

Using it from a Makefile:

```bash
export PKG_CONFIG_PATH=/home/j2c/Development/DFTB/dftbplus-25.1/install-qtes/lib/pkgconfig
pkg-config --cflags dftbplus        # -I.../modfiles plus conda flags
pkg-config --libs --static dftbplus # -ldftbplus plus OpenBLAS and runtime libs
```

The module files only work with the same compiler version (gfortran 15.2
from `rmg-dev`).

## Check (step 1)

A CH4 SCC single point with `bin/dftb+` (3ob-3-1 C/H parameters from
`/home/j2c/Development/DFTB/SlaKo/3ob-3-1/`, Fermi 6000 K, SCC tolerance
1e-6), run in `tests/runs/dftbp-step1-ch4_*`:

- converged in 4 SCC iterations; Mermin free energy −3.2291223287 Ha;
- forces with the expected tetrahedral symmetry (C: 0; each H:
  ±0.00267415 Ha/bohr per component).

## Notes

- The program banner reads "development version (commit: )" because the
  exported source has no git metadata; the library reports version 25.1
  (`dftbplus.pc`) and API 0.4.0.
- `dftb+` warns if the stack size is limited; use `ulimit -s unlimited`, as
  for `qtes.x`.
- By default `dftb+` writes `detailed.out`, `band.out`, `charges.bin` and
  `dftb_pin.hsd`. The driver must switch these off (audit §8.1).
- Forces: with ParserVersion 14 the keyword is `Analysis { PrintForces = Yes }`.
  It also makes the API compute forces when there is no `Driver` (see
  `README.md`).
- The conda environment adds its own optimization flags (e.g.
  `-march=nocona -mtune=haswell`) to the DFTB+ build through its
  `FFLAGS`/`CFLAGS`.
