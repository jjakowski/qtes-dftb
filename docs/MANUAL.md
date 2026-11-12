# qtes-dftb User Manual

**Status:** draft, step 1 (2026-11-05). Chapters marked **[to be written]**
need input from the code authors, mainly the physics background. All other
statements were checked against the source code and the example runs.

---

## Contents

1. [Introduction](#1-introduction)
2. [Theory in brief](#2-theory-in-brief) **[to be written]**
3. [Installation](#3-installation)
4. [Quick start](#4-quick-start)
5. [Input reference](#5-input-reference)
6. [Run modes](#6-run-modes)
7. [Running](#7-running)
8. [Output reference](#8-output-reference)
9. [Examples](#9-examples)
10. [Verification and regression tests](#10-verification-and-regression-tests)
11. [Troubleshooting](#11-troubleshooting)
12. [Developer notes](#12-developer-notes) **[to be written]**

Appendices: [A. Units](#a-units) · [B. Glossary](#b-glossary) ·
[C. Compile-time limits](#c-compile-time-limits)

---

## 1. Introduction

`qtes-dftb` simulates the dynamics of a molecular system in which a few
light nuclei (typically one H atom) are treated **quantum mechanically**
and all other atoms **classically**. The quantum nuclei are represented by
an ensemble of **quantum trajectories** (QTs) that sample their
wavepacket. Their quantum effects enter through an approximate quantum
potential, computed with the **linearized quantum force** (LQF) method. The
forces on all atoms come from **DFTB** electronic-structure calculations
done on the fly for every trajectory. The classical atoms move either in
the trajectory-averaged (Ehrenfest) force or independently in each
trajectory.

The program is a parallel Fortran code (`bohmd`, version 17/08/2016). It
distributes the trajectories over MPI ranks.

**Typical application:** an H atom scattering from, or chemisorbed on, a
hydrocarbon or a graphene flake. All included examples are of this type.

**Please cite:**

- S. Garashchuk, J. Jakowski, L. Wang, B. G. Sumpter, "Quantum
  Trajectory–Electronic Structure Approach for Exploring Nuclear Effects in
  the Dynamics of Nanomaterials", *J. Chem. Theory Comput.* **9**,
  5221–5235 (2013). doi:10.1021/ct4006147 (method)
- L. Wang, J. Jakowski, S. Garashchuk, "Adsorption of a Hydrogen Atom on a
  Graphene Flake Examined with Quantum Trajectory/Electronic Structure
  Dynamics", *J. Phys. Chem. C* **118**, 16175–16187 (2014).
  doi:10.1021/jp503261k (H on graphene application)

Related: L. Wang, J. Jakowski, S. Garashchuk, B. G. Sumpter, *J. Chem.
Theory Comput.* **12**, 4487–4500 (2016), doi:10.1021/acs.jctc.6b00126
(nonlinear quantum corrections, not implemented in this code version).

---

## 2. Theory in brief

**[to be written, with the code authors]** Planned content:

- the quantum-trajectory (Bohmian) formulation: QT positions, momenta,
  action `s` and weights `w`;
- the linearized quantum force (LQF) and the fit used in `qpot`;
- the Ehrenfest treatment of the classical atoms;
- the initial wavepacket and its sampling;
- the definitions of the energies and of the correlation function written
  to `ener` and `cor`.

These should follow the notation of the 2013 *JCTC* paper.

Facts already established from the code, to be built on:

- **Wavepacket:** the initial wavepacket of each quantum coordinate is a
  Gaussian with width parameter `a0`, centre `q0` and momentum `p0`. The
  sampling spread is `σ = 1/sqrt(4 a0)`.
- **Quantum-force fit:** the basis is linear, `(x, y, z, 1)` per quantum
  atom. The fit solves the weighted moment matrix of the QT positions.
- **Propagation:** velocity Verlet. Each large step `dt` is split into
  `kout` substeps.

---

## 3. Installation

### 3.1 Requirements

| Item | Tested with | Notes |
|---|---|---|
| Fortran compiler + MPI wrapper `mpif90` | gfortran 15.2 + OpenMPI (conda env `rmg-dev`) | Intel `ifort` + MKL also worked historically (settings kept as comments) |
| BLAS + LAPACK | OpenBLAS 0.3.30 | provides both; linked as `-lopenblas` |
| OpenMP | gfortran `-fopenmp` | |

On the development machine, everything comes from the `rmg-dev` conda
environment. Activate it before you build *and* before you run, because
`qtes.x` loads OpenBLAS and the MPI libraries from that environment at run
time.

### 3.2 Build

```bash
(cd dftb && make libdftb)   # 1. DFTB library  -> dftb/libdftb.a
make                        # 2. main program  -> qtes.x
```

The top-level `make` does not build the library, so step 1 is required
first, and again after any change in `dftb/`.

**Compiler flags:**

- `-O2 -fopenmp -fallow-argument-mismatch -std=legacy`, set in `Makefile`
  (`switch`, `LIB`) and in `dftb/GNUmake.inc` (`FOPT`).
- `-fallow-argument-mismatch` is required with gfortran 10 and later. The
  legacy code passes arrays of different types and ranks to the same
  routines, which newer gfortran rejects as errors unless this flag is set.
  With it, they become warnings, and the many warnings during the build are
  expected.
- To use Intel/MKL instead, restore the commented lines `-mkl=sequential`
  and `-openmp` in both files.

**Clean rebuild:**

```bash
make allclean && (cd dftb && make clean && rm -f *.o libdftb.a)
```

`make allclean` also deletes `fort.*` and `*.log` in the top directory.

### 3.3 Debug build

To locate a floating-point error, build a separate copy with:

```
-O0 -g -fbacktrace -ffpe-trap=zero,invalid -fallow-argument-mismatch -std=legacy -fopenmp
```

The program then stops at the first division by zero or invalid operation
and prints a backtrace with file and line numbers.

---

## 4. Quick start

```bash
conda activate rmg-dev                     # environment with mpif90/OpenBLAS
(cd dftb && make libdftb) && make          # build
tests/run_example.sh h-c4h9-14 4           # ~1 s smoke test on 4 MPI ranks
tests/run_example.sh h-c4h9-14-ref 4       # main regression test
R=$(ls -d tests/runs/h-c4h9-14-ref_* | tail -1)
diff -q $R/ener examples/h-c4h9-14-ref/reference/ener && echo "matches reference"
```

`tests/run_example.sh <example> [ranks]` does the following:

- creates a fresh directory `tests/runs/<example>_<timestamp>/`;
- copies in the example's inputs, the Slater–Koster files, any restart
  checkpoint, and `qtes.x`;
- sets `OMP_NUM_THREADS=1` (unless already set) and
  `OPENBLAS_NUM_THREADS=1`, and raises the stack limit;
- runs `mpirun -np <ranks> ./qtes.x > stdout.log`.

Never run `qtes.x` inside `examples/`. The program writes its output into
the current directory and would overwrite the reference data.

---

## 5. Input reference

The program takes **no command-line arguments**. It reads these files from
the current directory, under these exact names:

| File | Content | Read by |
|---|---|---|
| `INdyn` | dynamics settings + one line per atom (wavepacket, displacement, momentum, mass) | `read_input_mpi.f` |
| `INxyz` | model and analysis parameters + reference geometry | `read_input_mpi.f` |
| `IN_dftb` | DFTB atom types and Slater–Koster file names | `dftb/dftb_init_mpi.f` |
| `*.spl` | Slater–Koster parameter tables | DFTB library |
| `traj.sav` (+ old outputs) | only for a restart | `read_input_mpi.f`, `set_traj.f` |

All files are read with Fortran list-directed input. Values are separated
by spaces, and text after the last value on a line, such as the `\comment`
used in the examples, is ignored. Lines 1–9 of `INdyn` must each hold
exactly the values listed below, except where values are marked optional.

### 5.1 `INdyn`

| Line | Values | Meaning | Units, notes |
|---|---|---|---|
| 1 | `nat` | total number of atoms | must equal the number of geometry lines in `INxyz` |
| 2 | `nq` | number of **quantum** atoms; they must be listed **first** | examples: 1 (the H) |
| 3 | `kmax` | number of large time steps | on restart, reduced automatically by the time already done |
| 4 | `kout` | substeps per large step | |
| 5 | `dt` | large time step | a.u. of time; the integration step is `dt/kout` |
| 6 | `ipot isr` | `ipot`: quantum potential (0 off, 1 3-D LQF, 2 LQF along z). `isr`: initial sampling (0 off, 1 Gaussian) | `ipot > 0` requires `isr = 1` (§6) |
| 7 | `pow` | sampling cutoff: keep a sampled point if its Gaussian weight exceeds `10^-pow` | used only with `isr = 1` |
| 8 | `ipf` | 0 = Ehrenfest (averaged classical forces); any other value = independent classical copies | §6 |
| 9 | `np ist [nqtprt [nclprt]]` | `np`: number of QTs. `ist`: 0 = fresh start, 1 = restart from `traj.sav`. Optional print counts, default 20 each, `≤ 0` = all | §6.5, §8 |
| 10+ | one line per atom, in the same order as `INxyz` (blank lines are skipped) | see below | |

**Atom lines:**

```
quantum atom:    sym  a0   x0 y0 z0   px py pz   m
classical atom:  sym  (ignored)  x0 y0 z0   (px py pz ignored)   m
```

| Field | Meaning | Units, notes |
|---|---|---|
| `sym` | element symbol (informational) | |
| `a0` | width parameter of the initial Gaussian, the same for x, y, z | bohr⁻²; sampling spread `σ = 1/sqrt(4 a0)` |
| `x0 y0 z0` | displacement added to the `INxyz` position | bohr |
| `px py pz` | initial momentum of the wavepacket centre | a.u. Read for every atom, but **used only for quantum atoms**; classical atoms always start at rest |
| `m` | mass | electron masses (H = 1836, C = 22032). The `h-c4h9-14*` examples use masses ×10⁸ for the classical atoms, which effectively freezes the target |

### 5.2 `INxyz`

| Line | Values | Used in DFTB runs? | Meaning |
|---|---|---|---|
| 1 | `ql qls` | no | lattice constant and interaction cutoff of the old analytic model potentials (`derivs0/1/3`) |
| 2 | `dw sc` | no | force constant and coupling (`sc` is multiplied by `dw`) of the model potentials |
| 3 | `dch zh ch` | **`ch` only** | `dch`, `zh`: model-potential parameters. `ch`: reference z (bohr) used for the bonding probability in `aver` (§8.2) |
| 4 | `esh` | **yes** | energy shift subtracted from every DFTB energy, `V = E_DFTB − esh` (hartree) |
| 5 | `bnd` | **yes** | threshold for the bonding probability in `aver` (bohr) |
| 6 | `qmn qmx` | **yes** | z-window for the quantum atom (bohr), §6.4 |
| 7+ | `sym x y z itype` | yes | reference geometry, one line per atom, **in Ångström** (converted to bohr on input with 0.52918 Å/bohr). `itype` is the DFTB atom type, an index into `IN_dftb`. Blank lines before the geometry are allowed; lines after the first `nat` atoms are ignored |

The program builds the actual starting positions from the `INxyz`
geometry plus the `INdyn` displacements. Make sure no two atoms end up on
top of each other: DFTB then fails at the first step with
`ewevge: ier = …`. See `KNOWN_ISSUES.md` F2 and §11.

### 5.3 `IN_dftb`

```
2                ! ntype: number of atom types
2 1              ! lmax per type: number of angular shells (1 = s, 2 = s+p); type 1 = C, type 2 = H
'cc.spl'         ! SK files for every ordered pair (1,1), (1,2), (2,1), (2,2)
'ch-g2ms.spl'
'hc-g2ms.spl'
'hh.spl'
```

- **Slater–Koster files:** `examples/sk/` contains `cc.spl`, `hh.spl` and
  three C–H parameterizations: `ch/hc-g2ms.spl` (used by all examples),
  `ch/hc-b3.spl`, and `ch/hc.spl`.
- **Hard-coded DFTB settings** in `derivs` (`qtes.f`):
  - electronic temperature `telec = 6000 K` (Mermin free energy);
  - neutral system (`qnel = 0`);
  - SCC tolerance `1e-6`;
  - `ntype = 2`.

  Systems with more than two element types therefore need a code change.

### 5.4 Restart files

With `ist = 1` on line 9:

- `traj.sav` supplies the time and all QT positions, momenta, actions and
  weights.
- The output files `traj aver qt.dat pfit qtx qty qtz ener clx cly clz wft
  cor` must already exist. They are opened with `status='old'` and
  appended to.
- Keep all other inputs, including the print options, the same as in the
  original run. See example `h-c4h9-14-restart`.

---

## 6. Run modes

This section is a summary. `docs/DYNAMICS.md` describes each switch in
detail, with the code locations.

| Mode | Switch | Typical use | Requirements |
|---|---|---|---|
| Ehrenfest | `ipf = 0` | classical target responds to the averaged quantum atom | default in all examples |
| Independent classical copies | `ipf ≠ 0` | each QT carries its own target | only differs from Ehrenfest when `isr = 1` |
| No quantum potential | `ipot = 0` | quantum atom moves classically (ensemble of classical trajectories) | — |
| 3-D LQF | `ipot = 1` | quantum effects (zero-point energy, spreading) in x, y, z | `isr = 1` and enough QTs (≥ 4 with `w > 0`, in practice many more) |
| LQF along z | `ipot = 2` | quantum effects along the surface normal only | `isr = 1` |
| Sampling | `isr = 1`, `pow` | wavepacket represented by randomly sampled QTs | — |
| z-window | `qmn`, `qmx` | below `qmn`: QT removed from the ensemble; above `qmx`: DFTB skipped (free flight) | choose `qmx` where the potential is flat |
| Restart | `ist = 1` | continue a run from `traj.sav` | §5.4 |
| Print options | `nqtprt`, `nclprt` | width of `qtx..clz` | §8 |

**[to be written]** Physical guidance on choosing modes and parameters,
for example how many QTs, `a0`, `dt`, and `qmn`/`qmx`, with the code
authors.

---

## 7. Running

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ulimit -s unlimited                  # large local arrays live on the stack
mpirun -np <ranks> ./qtes.x > stdout.log 2>&1
```

- **Parallelization:**
  - Every rank holds a full copy of all QTs.
  - The DFTB calculations are split into contiguous blocks of
    `ceil(np/ranks)` QTs per rank.
  - The forces are gathered and broadcast after every substep.
  - **Rank count rule:** **every rank must get at least one QT**, i.e.
    `(ranks-1)·ceil(np/ranks) < np`. Otherwise the run **hangs** at the
    start, with no error message (`KNOWN_ISSUES.md` #8).
  - Uneven splits are fine. For example, 48 QTs run on 4, 5, 7 or 10
    ranks with bit-identical results, but hang on 13 ranks.
  - The safest choice is a divisor of `np`. Results do not depend on the
    rank count.
- **Threads:** keep `OPENBLAS_NUM_THREADS=1`. This matches the original
  `-mkl=sequential` and avoids oversubscribing cores with MPI ranks ×
  BLAS threads.
- **Stack:** `ulimit -s unlimited` is needed. With `-fopenmp`, large local
  arrays (about 28 MB in `trot`) are allocated on the stack.
- **Cluster:** `docs/sbatch.script` is an old SLURM example (24 tasks per
  node, `OMP_NUM_THREADS=1`).

**Measured run times** (gfortran + OpenBLAS, this machine):

| Example | Atoms | QTs | Substeps | Ranks | Wall time |
|---|---|---|---|---|---|
| `h-c4h9-14` | 14 | 48 | 8 | 4 | ~1 s |
| `h-c4h9-14-ref` | 14 | 48 | 50 | 4 | ~3 s |
| `h-flake-53` | 53 | 48 | 30 | 48 | ~4 s |
| `h-flake-111` | 111 | 12 | 50 | 12 | ~25 s |

---

## 8. Output reference

All output files are written by MPI rank 0 into the current directory.

- **Units:** time in a.u., positions in bohr, energies in hartree.
- **Columns:** in all files except `wf0`, `wft`, `qt.dat` and `traj.sav`,
  column 1 is the time.
- **When written:**
  - **per substep:** `ener`, `cor`;
  - **per large step:** `aver`, `traj`, `qtx..clz`, `pfit`, plus `t = 0`
    for most of them;
  - `qt.dat`: per large step, without `t = 0`;
  - `traj.sav`: overwritten at every large step.

### 8.1 Files

| File | Content and columns |
|---|---|
| `ener` | `t, E_total, <Q>, <V>, <K_quantum>, <K_classical>`: weighted averages over the QTs; `V` includes the `esh` shift. Energy conservation of `E_total` is the main health check |
| `cor` | `2t, Re C, Im C`: correlation function `C = Σ w exp(2 i s) / Σ w` (an autocorrelation if ψ(0) is real) |
| `aver` | `t, norm, P, <x_1> … <x_40>`: norm = `Σ w` (drops below 1 when QTs leave below `qmn`); `P` = bonding probability (§8.2); weighted average of the first 40 coordinates |
| `traj` | `t`, the first `4·nqd` coordinates of QT 2, and `s(2)`: a representative trajectory |
| `qtx qty qtz` | `t`, then x/y/z of the quantum atom in QT 1 … `min(np, nqtprt)` |
| `clx cly clz` | `ipf = 0`: `t`, then x/y/z of classical atoms 1 … `min(ncl, nclprt)` from the first QT with `w > 0`. `ipf ≠ 0`: the first classical atom in QT 1 … `min(np, nqtprt)` |
| `pfit` | `t` and the LQF fit coefficients (`ipot > 0`; empty otherwise) |
| `qt.dat` | 2 lines per QT per large step: all positions + `w`, then all momenta + `s` (full precision) |
| `traj.sav` | restart checkpoint: time, then the same 2 lines per QT |
| `wf0`, `wft` | initial and final positions of every QT (first 40 coordinates). On restart, `wf0` holds the restart state and `wft` is appended to |
| `trajt`, `dftbt` | timings: start, end, elapsed. They differ on every run |
| `out.dftb` | echo of the DFTB setup read from `IN_dftb` |
| `stdout.log` | screen output: setup summary, including `PRINT OPTIONS`, initial energies, step counter |
| `fort.1000+rank` | per-rank trace written at initialization |

`docs/sg-notes-files` gives the column formats in the original notation.

### 8.2 Bonding probability `P` (column 3 of `aver`)

`P = Σ w_i [ |z_H,i − ch| ≥ bnd ] / Σ w_i`, using the **absolute z
coordinate** of the quantum atom and `ch`, `bnd` from `INxyz`.

With the values in the current examples (`ch = 0.65` bohr, `bnd = 1`
bohr), the H is always more than `bnd` away from `z = ch`, so `P = 1`
throughout. **[to be clarified with the authors]** What `P` is meant to
measure. A commented-out older line used the H–C distance instead of the
absolute z.

---

## 9. Examples

All examples are in `examples/`. `examples/README.md` describes each one
in detail, including how its reference data was made.

| Example | System | Purpose |
|---|---|---|
| `h-c4h9-14-ref` | H + C₄H₉ (frozen), 48 QTs | main regression test |
| `h-c4h9-14-restart` | same, restart from t = 80 | tests the restart path |
| `h-c4h9-14` | same, short | ~1 s smoke test |
| `h-flake-53` | H + 52-atom graphene flake | larger system, Ehrenfest target |
| `h-flake-53-print` | same, `nqtprt = 48`, `nclprt = 0` | demonstrates the print options |
| `h-flake-111` | H chemisorbed on a 110-atom flake, `ipot = 1`, `isr = 1` | quantum potential with sampling |

The graphene-flake examples correspond to the system studied in Wang,
Jakowski, Garashchuk, *J. Phys. Chem. C* **118**, 16175 (2014).

**[to be written]** Step-by-step tutorials with plots.

---

## 10. Verification and regression tests

- **Reproducibility:** runs are bit-reproducible, and the results are
  the same for any valid rank count. This was tested with 4, 5, 7 and 10
  ranks on `h-c4h9-14-ref`.
  Every example has a `reference/` folder.
- **After a code change:** run all examples. Every file except
  `trajt`/`dftbt` must be byte-identical to the reference, unless the
  change is meant to alter results.

```bash
for ex in "h-c4h9-14-ref 4" "h-c4h9-14-restart 4" "h-flake-53 48" \
          "h-flake-53-print 48" "h-flake-111 12"; do
  set -- $ex; tests/run_example.sh $1 $2 >/dev/null
  R=$(ls -d tests/runs/$1_2* | tail -1)
  for f in ener cor aver qt.dat traj.sav qtx clx; do
    cmp -s $R/$f examples/$1/reference/$f || echo "$1: $f differs"
  done
done
```

- **Physics checks:**
  - `E_total` in `ener` should be conserved, to all printed digits in
    the examples.
  - A restart from t = 80 reproduces a single run from 0 to 100 to
    round-off.
  - With per-rank checksums, all ranks stay bit-identical (see
    `KNOWN_ISSUES.md` F7).

---

## 11. Troubleshooting

| Message or symptom | Cause | Fix |
|---|---|---|
| `Bad integer for item 2 in list input` at `read_input_mpi.f` | a line of `INdyn` has too few values, e.g. line 6 needs `ipot isr` | give all required values (§5.1) |
| `ewevge: ier = …` in `fort.10x`, run stops before t = 0 | DFTB failed, usually because two atoms overlap after the `INdyn` displacements | check the starting distances (`INxyz` + displacements) |
| `QPOT ERROR … only N QTs with w>0; fit needs M` | too few QTs left inside the `qmn` window for the quantum-force fit | more QTs, a lower `qmn`, or `ipot = 0` |
| `QPOT ERROR … DPOSV FAILED … (isr=0?)` | the fit matrix is singular: all QTs identical (`isr = 0`) or lying in a plane | use `isr = 1` with `ipot > 0` |
| energies grow quickly before a `QPOT ERROR` | the fit is ill-conditioned with few, clustered QTs | more QTs; keep them inside the window |
| restart stops opening a file | restart needs the previous output files | copy them into the run directory (§5.4) |
| `IEEE_DENORMAL` note at exit | harmless: products of near-zero coordinates | none |
| `MPI_FINALIZE ERR 0` (one per rank) | the return code of `MPI_Finalize`; 0 means success | none |
| step in `E_total` when a QT re-enters below `qmx` | free flight with a frozen potential above `qmx` | raise `qmx` to where the potential is flat |
| run hangs right after `FINISHED SAMPLING`, all ranks at 100 % CPU | some MPI rank has no QTs (`KNOWN_ISSUES.md` #8) | use a rank count with every rank getting ≥ 1 QT, e.g. a divisor of `np` |
| run hangs or crashes at start with a large system | stack too small | `ulimit -s unlimited` |

Open and fixed bugs are tracked in `docs/KNOWN_ISSUES.md`.

---

## 12. Developer notes

**[to be written]** Planned content:

- the source layout (`qtes.f` subroutines `init_traj`, `trot`, `qpot`,
  `over`, `aver`, `derivs`; `read_input_mpi.f`; `set_traj.f`; `dftb/`);
- the replicated-data MPI design and its rules (all ranks must make
  identical decisions; see F7);
- the fixed-form Fortran 72-column limit;
- where to add a new potential (`derivs0..3`) or an output;
- the test procedure before committing.

---

## A. Units

| Quantity | Unit | Conversion |
|---|---|---|
| length | bohr (a₀) | 1 bohr = 0.52918 Å (value used in the code) |
| time | a.u. | 1 a.u. = 0.024189 fs |
| energy | hartree (Eh) | 1 Eh = 27.211 eV |
| mass | electron mass (mₑ) | H = 1836, C = 22032 (12 × 1836) |
| momentum | a.u. | ħ/a₀ |

## B. Glossary

| Term | Meaning |
|---|---|
| QT | quantum trajectory: one member of the ensemble that samples the quantum nuclei |
| QDOF / `nqd` | quantum degrees of freedom: 3 × number of quantum atoms |
| LQF | linearized quantum force: the quantum potential from a linear fit over the QTs |
| Ehrenfest | the classical atoms move in the QT-averaged force (`ipf = 0`) |
| `w`, `s` | weight and action (phase) of a QT |
| z-window | `qmn`/`qmx`: the z-range of the quantum atom inside which DFTB is computed |

## C. Compile-time limits

| Parameter | File | Value | Limits |
|---|---|---|---|
| `ndx` | `bohm.inc` | 700 | total coordinates, i.e. at most 233 atoms |
| `npx` | `bohm.inc` | 5000 | QTs |
| `nqx` | `bohm.inc` | 200 | quantum coordinates |
| `nndim` | `dftb/maxima.inc` | 800 | atoms in DFTB |
| `mdim` | `dftb/maxima.inc` | 3200 | DFTB orbitals |
| `maxtyp` | `dftb/maxima.inc` | 6 | DFTB atom types (`derivs` currently passes 2) |

After changing a limit, do a clean rebuild of both `dftb/` and the main
program.
