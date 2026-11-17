# Audit: replacing the built-in DFTB code with DFTB+

**Date:** 2026-11-17. **Status:** audit done, decisions recorded (§11). No
code was changed, and nothing was built or run in either repository.

**Question:** can `qtes-dftb` call DFTB+, maintained externally, as a library
to get the energy and forces, instead of the old DFTB code in `dftb/`? What
would it take?

**First milestone, as defined:** a small Fortran driver, in the style of
`qtes-dftb`, that passes coordinates and input parameters to DFTB+ and gets
the energy and forces back. Reproducing results of the old code is not
required. Physical settings (electronic temperature, SCC tolerance, ...) must
be user input, not hard-coded.

---

## 1. Summary

**Yes, this is possible, and the first milestone is low risk.**

- DFTB+ ships a supported Fortran library API (module `dftbplus`, type
  `TDftbPlus`), plus a C API. One of its own test programs
  (`test/src/dftbp/api/mm/testers/test_treeinit.f90`) does exactly what
  milestone 1 needs:
  1. create an instance;
  2. give it the input;
  3. call `setGeometry(coords)`;
  4. call `getEnergy(merminEnergy)` and `getGradients(gradients)`;
  5. repeat steps 3–4 for a new geometry.
- **Units and quantities match `qtes-dftb`:**
  - coordinates go in, in bohr;
  - the energy comes back in hartree and is the Mermin free energy, the same
    quantity the old code returns (`Emermin`);
  - the gradients come back in hartree/bohr, which is what `qtes-dftb` uses as
    `dv`.
- **All physical settings come from a DFTB+ input** in its usual HSD format.
  This can be a file, or built in code, or both. Nothing has to be hard-coded,
  which answers the request to make `telec`, the SCC tolerance, etc.
  user-adjustable.
- **The `qtes-dftb` code has a single, clean connection point.** All DFTB work
  goes through `derivs(nat, iat, C, Energy, frc)` in `qtes.f`, and the old
  library is set up once by `dftb_init_mpi`. An adapter behind `derivs` can
  switch engines without touching the dynamics.
- **Building:** a minimal DFTB+ library build needs only CMake, a Fortran
  compiler, Python (for its bundled preprocessor) and BLAS/LAPACK. All of these
  are already in the `rmg-dev` environment (CMake 4.1.2, gfortran 15.2, Python
  3.12, OpenBLAS). The minimal build needs no extra downloads.
- **The main things to handle:**
  - turn off the files DFTB+ writes on every call;
  - SCC charges carry over between calls on the same rank;
  - one DFTB+ instance per MPI process, unless an instance-safe build is used;
  - DFTB+ stops the program on errors.

  None of these blocks milestone 1 (details in §8).

---

## 2. What was examined

| Item | Value |
|---|---|
| DFTB+ checkout | `/home/j2c/Development/DFTB/dftbplus`, `25.1-17-g7e1a071d` (commit `7e1a071d`), clean working tree |
| API version | 0.4.0 (`src/dftbp/api/mm/API_VERSION`) |
| Licence | LGPL-3 (`COPYING.LESSER`). Linking from `qtes-dftb` is allowed; changes to DFTB+ itself would have to stay LGPL |
| API sources | `src/dftbp/api/mm/`: `mmapi.F90` (Fortran API, 1335 lines), `capi.F90` and `dftbplus.h` (C API) |
| API usage examples | `test/src/dftbp/api/mm/testers/` (28 test programs, Fortran and C) |
| Documentation | `INSTALL.rst` ("Using DFTB+ as a library"), `doc/dftb+/manual/api.tex`, `config.cmake` |
| `qtes-dftb` side | `qtes.f` (`derivs`), `dftb/do_dftb.f`, `dftb/dftb_init_mpi.f`, `examples/*/IN_dftb`, `examples/sk/*.spl` |

The machine also has other DFTB+ checkouts with existing library builds:
`../dftbplus_tblite/_build`, and `../dftbplus_forces_print/my_build` and
`build_tblite_sdftd3_gfortran`. Each provides `libdftbplus.a`, `dftbplus.pc`
and a CMake package config. They were not examined further, and the plan
below does not depend on them.

---

## 3. The current DFTB interface in `qtes-dftb`

### 3.1 Call chain

```
trot / init_traj / set_traj            (one call per QT, on its own rank)
  └─ derivs(nat, iat, C, Energy, frc)  qtes.f
       └─ do_dftb(...)                 dftb/do_dftb.f (old DYLCAO code, 1991)
setup, once:  dftb_init_mpi(ntype)     dftb/dftb_init_mpi.f
```

### 3.2 What `derivs` passes and gets back

| Item | Meaning | Units |
|---|---|---|
| `nat` | number of atoms | |
| `iat(nat)` | atom type per atom (index into `IN_dftb`: 1 = C, 2 = H in the examples) | |
| `C(3,nat)` | coordinates | bohr |
| `Energy` (out) | Mermin free energy `E - T·S` (`Emermin`) | hartree |
| `frc(3,nat)` (out) | `do_dftb` returns forces (`-gradient`), and `derivs` flips the sign, so **`derivs` returns the gradient** | hartree/bohr |

After the call, `qtes-dftb` subtracts `esh` (`INxyz` line 4) from the
energy.

### 3.3 Settings hard-coded in `derivs` (to become input)

| Setting | Value now | DFTB+ input equivalent |
|---|---|---|
| electronic temperature | `telec = 6000.0` K | `Hamiltonian = DFTB { Filling = Fermi { Temperature [K] = 6000 } }` |
| total charge | `qnel = 0` (neutral) | `Charge = 0.0` |
| SCC on/off | `scfhelp = .true.` | `Scc = Yes` |
| SCC tolerance | `scftol = 1e-6` | `SccTolerance = 1e-6` |
| number of types | `ntype = 2` | `TypeNames` in `Geometry` (any number) |
| angular momentum per type | `lmax` from `IN_dftb` line 2 (`2 1` = C up to p, H s only) | `MaxAngularMomentum { C = "p"; H = "s" }` |
| SK files per type pair | `IN_dftb` lines 3+ | `SlaterKosterFiles { C-C = "cc.spl"; C-H = "ch-g2ms.spl"; ... }` |

### 3.4 State between calls

- The Slater–Koster tables are read once: rank 0 reads them in
  `dftb_init_mpi`, and they are broadcast through a buffer.
- **The SCC charges are not kept.** Each call starts from neutral atoms
  (`qmat(i) = qzero(izp(i))`, `do_dftb.f:250`). So every call is independent
  of the previous one, which makes results independent of which rank computes
  which QT.

### 3.5 Minor observation (not a bug to fix now)

- `derivs` passes an undefined scalar `V`, and leaves out the last argument
  `IV`, in its call to `do_dftb(..., telec, V, IV)`.
- `do_dftb` only declares `V` and `IV` and never uses them, so this is
  harmless. It is one of the argument mismatches that need
  `-fallow-argument-mismatch`.
- It would disappear with the old engine.

---

## 4. The DFTB+ library API

### 4.1 Calls milestone 1 needs (Fortran, `use dftbplus`)

| Call | Purpose |
|---|---|
| `TDftbPlus_init(dftbp [, outputUnit, mpiComm, devNull])` | create an instance. `outputUnit` sends DFTB+ screen output (and errors) to a unit of our choice. `mpiComm` only for MPI builds of DFTB+ |
| `dftbp%getInputFromFile("dftb_in.hsd", input)` | read a normal DFTB+ input file into an input tree |
| `dftbp%getEmptyInput(input)` + `setChild` / `setChildValue` | build or extend the input tree in code (`test_treeinit.f90`) |
| `dftbp%setupCalculator(input)` | process the input, read the SK files, allocate |
| `dftbp%setGeometry(coords [, latVecs])` | new coordinates, `coords(3,nAtom)` in bohr. The atom order and species stay fixed |
| `dftbp%getEnergy(merminEnergy)` | runs the calculation if needed; returns the Mermin free energy (hartree) |
| `dftbp%getGradients(gradients)` | `gradients(3,nAtom)`, hartree/bohr |
| `dftbp%getGrossCharges(q)` | Mulliken charges, optional diagnostics |
| `TDftbPlus_destruct(dftbp)` | release the instance |

A typical sequence for `qtes-dftb`: create the instance and set up the
calculator **once per MPI rank**. Then, for every QT the rank owns and every
substep, call `setGeometry`, then `getEnergy` and `getGradients`.

### 4.2 Input: file, code, or both

- **From a file:** a normal `dftb_in.hsd` can hold all settings, including
  `Geometry`. This is the easiest route for milestone 1, and it lets
  users set any DFTB+ option without changes to `qtes-dftb`.
- **From code:** the input tree can be built entirely in code, as
  `test_treeinit.f90` does. It can also be read from a file and then
  extended. The driver could therefore add the `Geometry` block
  (`TypeNames`, species, start coordinates from `INxyz`) itself, and leave
  only the physics settings to the file. That avoids repeating the geometry
  in two input files.

### 4.3 Instances and parallelism

- **Default build (`INSTANCE_SAFE_BUILD = FALSE`):** one DFTB+ instance per
  process. A second `TDftbPlus_init` in the same process stops with an
  error. This fits `qtes-dftb`, where each rank needs one instance, reused
  for all its QTs.
- **Instance-safe build:** allows several instances per process. That would
  only be needed for one instance per QT, for example to keep each QT's
  charges (§8.2). It excludes some optional components.
- **MPI:**
  - A **serial** DFTB+ build (`WITH_MPI = FALSE`, the default) can be linked
    into an MPI program. Each rank then runs its own independent DFTB+,
    which is exactly the current model: every rank computes whole QTs.
  - An MPI build of DFTB+ would only matter for splitting one large DFTB
    calculation over several ranks. That is not needed for milestone 1.
- **OpenMP:** `WITH_OMP` is on by default. With many MPI ranks per node,
  either build without OpenMP or run with `OMP_NUM_THREADS=1`. This matches
  the current setup.

### 4.4 C API

`dftbplus.h` / `capi.F90` offer the same core calls from C. The Fortran
API is the natural choice here: `qtes-dftb` is Fortran, and both would be
compiled with the same gfortran. The C API would only be needed to avoid
depending on compiler-specific Fortran module files (§7.3).

---

## 5. Mapping old engine to DFTB+

| Need of `qtes-dftb` | Old engine | DFTB+ | Status |
|---|---|---|---|
| coordinates in bohr, `(3,nat)` | yes | `setGeometry(coords)`, bohr, `(3,nAtom)` | direct |
| Mermin free energy | `Emermin` | `getEnergy` returns the Mermin free energy | direct |
| gradient | `-Frc`, sign flipped in `derivs` | `getGradients` | direct (no sign flip) |
| electronic temperature | hard-coded 6000 K | `Filling = Fermi { Temperature [K] = ... }` | input |
| SCC, tolerance, charge | hard-coded | `Scc`, `SccTolerance`, `Charge` | input |
| element types | `iat` = index into `IN_dftb`, max 2 in `derivs` | `TypeNames` + species index per atom, any number | direct, limit removed |
| SK parameters | `.spl` files | `.skf` files (same format, §6) | to verify in the first test |
| neutral start each call | yes | no: the previous charges are the start guess (§8.2) | different, acceptable |
| `esh` energy shift | applied in `qtes-dftb` | unchanged | unaffected |
| reading parameters once | rank 0, then broadcast | each rank reads its own | fine (small files) |

---

## 6. Slater–Koster parameters

- **Same format.** The `.spl` files in `examples/sk/` use the classic SKF
  layout:
  1. grid spacing and number of points;
  2. on-site energies and Hubbard U (homonuclear files);
  3. mass and polynomial line;
  4. the H/S table;
  5. a `Spline` block for the repulsion.

  This agrees with your statement that `.spl` and `.skf` are compatible.
  DFTB+ takes any file name per type pair, so no renaming is needed:
  `C-H = "ch-g2ms.spl"` works.
- **Two things to check in the first test**, by comparing the DFTB+ result
  with the DFTB+ standalone binary on the same input:
  - The files use Fortran repeat counts (`19*0.0`, `20*1.0`).
  - The heteronuclear `ch-g2ms.spl` has two `20*1.0` lines at the top.

  If DFTB+ rejects either, a one-off conversion script fixes it.
- **No standard sets downloaded.** `external/slakos` in the DFTB+ checkout
  holds only a status file. They are not needed if the `.spl` files work.

---

## 7. Building and linking

### 7.1 Minimal DFTB+ library build

- Build **outside both repositories**, e.g. a build directory and an
  install prefix under `~/Development/DFTB/`.
- Build with the same compiler as `qtes-dftb` (gfortran 15.2 from `rmg-dev`).

```bash
cmake -S /home/j2c/Development/DFTB/dftbplus -B <build-dir> \
      -DWITH_API=ON -DWITH_MPI=OFF -DWITH_OMP=OFF \
      -DCMAKE_INSTALL_PREFIX=<install-dir> \
      -DLAPACK_LIBRARY=openblas -DBLAS_LIBRARY=openblas
cmake --build <build-dir> -j && cmake --install <build-dir>
```

(The exact flags are to be confirmed when building; `WITH_API` is already
on by default.)

**What the minimal build needs:**
- **Required:** CMake ≥ 3.16, Fortran and C compilers, Python (to run the
  bundled `fypp` preprocessor), BLAS/LAPACK.
- **Bundled:** `fypp` and `xmlf90` are in `external/`.
- **Not needed for this build:** `toml-f` (only with
  `INCLUDE_INDIRECT_DEPS`), `mpifx` and `scalapackfx` (only with MPI). The
  uninitialised git submodules in `external/` all belong to optional
  components (tblite, dftd4, s-dftd3, mbd, libnegf, chimes, ...). So no
  network access should be needed.

### 7.2 What the install provides

- `lib/libdftbplus.a` (or `.so` with `BUILD_SHARED_LIBS=ON`), plus the
  helper libraries it needs;
- the Fortran module files in `include/dftbplus/modfiles/`;
- the C header `dftbplus.h`;
- a pkg-config file `dftbplus.pc`, and a CMake package config
  (`dftbplus-config.cmake`).

The `qtes-dftb` `Makefile` can take the compile and link flags from
`pkg-config --cflags --libs dftbplus`, so no paths are hard-coded.

### 7.3 Compiler compatibility

- **Fortran module files** (`.mod`) are specific to the compiler and its
  version. DFTB+ and the driver must be compiled with the same gfortran.
  Building DFTB+ in `rmg-dev` guarantees this.
- **A pre-built DFTB+** (e.g. the conda-forge package) could work if its
  gfortran matches, but that is less controlled.
- **Flags:** the legacy `qtes-dftb` flags (`-std=legacy`,
  `-fallow-argument-mismatch`) only apply to the old fixed-form files. The
  new adapter would be a separate free-form `.f90` file compiled with normal
  flags. The fixed-form code just calls a plain subroutine in it.

### 7.4 Fortran 77 (`qtes-dftb`) and Fortran 90/2003 (DFTB+) together

This is not a problem in itself. gfortran compiles both into one program, and
code in one style can call code in the other. "Fixed form" is a source layout,
not a language level: gfortran accepts modern features (e.g. `use`) in
fixed-form files. Practical rules:

- **All modern Fortran stays in a small free-form adapter** (`.f90`). DFTB+
  types and objects (`TDftbPlus`, input trees), its module and its
  assumed-shape arrays live there. The old code only calls plain subroutines
  with plain arguments, e.g. `dftbp_energy_gradient(nat, C, E, grad)` with
  `C(3,nat)`, the same kind of call it makes to `do_dftb` today.
- **Argument checking.** Old-style calls are not checked by the compiler, so
  a wrong argument type would go unnoticed (implicit typing caused several
  earlier bugs). The adapter should have exactly typed, explicit-shape
  arguments. `derivs` can also `use` the adapter module, so that the
  compiler checks each call.
- **Number kinds must match.** `real*8` in `qtes-dftb` is DFTB+'s
  `real(dp)`, and both use default 4-byte integers. DFTB+ must not be built
  with `-fdefault-integer-8`, `-fdefault-real-8` or similar flags.
- **Separate compiler flags.** The legacy flags (`-std=legacy`,
  `-fallow-argument-mismatch`) apply only to the old files; the adapter gets
  normal flags. Everything is compiled and linked by `mpif90` (gfortran), so
  there is one Fortran runtime.
- **Same gfortran version** for DFTB+ and the adapter, because of the module
  files (§7.3).
- **Symbol clashes are unlikely.** DFTB+ routines live in modules, so their
  linker names carry a module prefix and cannot clash with old global names
  such as `fermi` or `broyden` in `libdftb.a`. This matters only if both
  engines are linked into one executable; to be confirmed in the first link
  test.

---

## 8. Points to watch

### 8.1 Files written by DFTB+ on every call

- By default DFTB+ writes `detailed.out`, `band.out` and `charges.bin` to
  the current directory (the defaults `WriteDetailedOut = Yes` and
  `WriteCharges = Yes` are in `dftbplus/parser.F90`). With thousands of
  calls, and all ranks in the same directory, that means overhead and
  ranks overwriting each other's files.
- **Fix, in the input:**
  ```
  Options { WriteDetailedOut = No; WriteCharges = No }
  Analysis { WriteBandOut = No }
  ```
- Also pass `outputUnit`, to send the DFTB+ screen output to one file per
  rank or to the null device. Which files are actually written in API mode
  should be confirmed in the first test.

### 8.2 SCC charges carry over between calls

- `setGeometry` marks the coordinates as changed but keeps the charges
  (`mainapi.F90`). So the next SCC cycle starts from the previous call's
  converged charges. On one rank, the previous call is usually a
  *different* QT (or the same QT one substep earlier).
- **Effect:** the converged result is the same to within `SccTolerance`,
  and convergence is usually faster than from neutral atoms. But results
  are no longer bit-identical between different rank counts, because they
  depend slightly on the order of calls.
- **This matters for testing:** the current regression tests rely on that
  bit-identity. With DFTB+, comparisons between rank counts need a
  tolerance, unless the start charges are reset before each call.
- **No files are involved.** The carried-over charges live in memory, inside
  the instance. Writing `charges.bin` is a separate default, switched off
  with `WriteCharges = No` (§8.1). Separate run folders per QT are therefore
  not needed. They would also not scale: 1000+ QTs writing and reading files
  in real time from many ranks would become a bottleneck.
- **The three possible start guesses:**

  | Option | How | DFTB+ 25.1 |
  |---|---|---|
  | A. the previous call's charges (default) | nothing to do | works |
  | B. each QT's own charges from its previous step, kept in memory | `qtes-dftb` stores one charge array per QT (1000 QTs × 233 atoms ≈ 2 MB) and passes it to DFTB+ before each call | reading charges works (`getGrossCharges`); **setting the SCC start charges is not in the API** |
  | C. neutral atoms every call, like the old code | reset the start guess before each call | **not in the API** (recreating the instance per call would be far too slow) |

  `setRefCharges` sets the *reference* (neutral-atom) charges, not the SCC
  start guess.
- **Decision** (§11): milestone 1 uses option A. Options B and C need an
  API call to set the SCC start charges (§8.3).

### 8.3 Plan: adding an API call to set the SCC start charges

We would develop the call ourselves and offer it upstream as a pull request,
following the workflow used for the `forces-print` feature (fork
`jjakowski/dftbplus`, a feature branch in its own worktree under
`~/Development/DFTB/`, notes in `~/Development/DFTB/readme-short.md`). This
gives control over the design and lets `qtes-dftb` use the call before it is
merged.

**Scope** (from the audit, to be confirmed in the code):
- **API layer:**
  - a new call in the Fortran API (`src/dftbp/api/mm/mmapi.F90`), e.g.
    `setSccCharges(atomCharges)`, the counterpart of `getGrossCharges`;
  - its implementation in `src/dftbp/dftbplus/mainapi.F90`;
  - the C wrapper in `capi.F90` and `dftbplus.h`.
- **Inside DFTB+:** the charges are stored per orbital and spin, and the
  SCC mixer (e.g. Broyden) keeps a history that must be reset when the start
  charges change. A first version could accept atomic charges only, for
  non-spin-polarized SCC, and stop with a clear error otherwise. Passing the
  neutral reference charges through the same call would also give option C.
- **Around the code:** an API test program in
  `test/src/dftbp/api/mm/testers/`, an API version bump (now 0.4.0), a
  `CHANGELOG.rst` entry, and a section in `doc/dftb+/manual/api.tex`.

**Approach:**
1. **Open a GitHub issue first**, describing the use case (one charge vector
   per trajectory kept in memory, thousands of calls per run), so that the
   DFTB+ developers can say how they would like it designed. This should
   save review rounds later.
2. Implement and test on the feature branch; build `qtes-dftb` against that
   install for testing.
3. Open the pull request.

**`qtes-dftb` keeps working with stock DFTB+ 25.1** (option A), so it never
depends only on an unmerged fork. The adapter uses the new call only when it
is available.

### 8.4 Errors stop the program

DFTB+ reports input or convergence errors through its `error` routine,
which stops the process. In an MPI run, one rank stopping can leave the
others hanging, like the issues fixed in F6 and F9. The adapter should
check its inputs before calling DFTB+. Error handling in API mode should
be tested, and `mpirun` should end the job when one rank exits.

### 8.5 Per-call cost

A run makes many calls: QTs × substeps × steps, spread over the ranks. The
API overhead per call (`setGeometry`, neighbour list update, SCC from the
previous charges) is unknown. It should be measured in milestone 1 against
the old engine on the same geometries.

### 8.6 Other

- `qtes-dftb` uses `ndx = 700` coordinates (about 233 atoms) and the old
  engine's limits (`dftb/maxima.inc`). DFTB+ allocates dynamically, so its
  own limits disappear, but `bohm.inc` still applies.
- **Choice of input file name:** DFTB+ reads `dftb_in.hsd` by default, and
  `getInputFromFile` takes any name. A name like `IN_dftbplus.hsd` would
  follow the `IN*` style of `qtes-dftb`.

---

## 9. Proposed design (for later steps)

1. **New input file** `IN_dftbplus.hsd`, in normal DFTB+ format, holding all
   DFTB settings (temperature, SCC, tolerance, charge, SK files, angular
   momenta, output switches). Nothing physical stays hard-coded. The
   geometry comes from `INxyz` as now: the adapter adds the `Geometry` block
   in code (`TypeNames` from `IN_dftbplus.hsd` or `INxyz`, species from the
   `iat` column).
2. **Adapter module** `dftbplus_engine.f90` (free-form, normal flags) with
   two plain subroutines callable from fixed-form code:
   - `dftbp_setup(nat, iat, C)`: once per rank;
   - `dftbp_energy_gradient(nat, C, E, grad)`: per call.
3. **Engine switch:** `derivs` calls either the old `do_dftb` or the
   adapter. The switch could be chosen by the presence of
   `IN_dftbplus.hsd`, or by an explicit input flag. Keeping the old engine
   at first allows side-by-side comparison on the same inputs.
4. **Build:** an optional `Makefile` target (e.g. `make DFTBPLUS=1`) that
   adds the adapter and the `pkg-config dftbplus` flags. The default build
   stays as it is.

---

## 10. Plan for milestone 1 (proposed, not started)

The current roadmap (milestones 1–6) and the status of each step are kept in
`agent_state.md` in the repository root.

| Step | Work | Where | Success check |
|---|---|---|---|
| 1 | Build and install the DFTB+ library (serial, API on, OpenBLAS) with `rmg-dev` gfortran | outside both repos, e.g. `~/Development/DFTB/build-qtes/` and `install-qtes/` | `dftbplus.pc` and module files installed; the standalone `dftb+` binary from the same build runs a test input |
| 2 | Standalone Fortran driver `dftbp_driver.f90` (about 100 lines): read a geometry (bohr) and `IN_dftbplus.hsd`, call `TDftbPlus_init` → `getInputFromFile` → add `Geometry` → `setupCalculator` → `setGeometry` → `getEnergy` / `getGradients`, print E and the gradients; then move one atom and repeat | new folder in `qtes-dftb`, e.g. `tools/dftbplus/`, with its own small Makefile | builds and links against the installed library |
| 3 | Run it on the `h-c4h9-14` geometry with the `.spl` files, `Temperature = 6000 K`, SCC on, tolerance 1e-6 | `tests/runs/` | (a) same energy and gradient as the `dftb+` binary on the same input; (b) gradient agrees with a finite-difference derivative of the energy; (c) same order of magnitude as the old engine (no exact match expected) |
| 4 | Several geometries in a row on one instance (a mini-trajectory, several QTs as on one rank), charge option A; measure per call the SCC iterations and the time, against the old engine; check how much E and forces depend on the order of calls | `tests/runs/` | stable SCC; cost per call known; dependence on call order within the SCC tolerance; no files written by DFTB+ |
| 5 | Optional: MPI version of the driver, one instance per rank, each rank computing its own geometries | `tests/runs/` | runs on several ranks with no file clashes |

After milestone 1 comes the adapter and the engine switch inside `qtes.f`
(§9), then a dynamics run with DFTB+ forces, checked for energy
conservation. That run would need new reference outputs.

**Rough effort:**
- milestone 1, steps 1–4: about one to two working days, most of it
  building and checking the SK files and the output settings;
- the integration into `qtes-dftb`, with tests and docs: a few more days.

---

## 11. Decisions and open questions

**Decided (2026-11-17):**

1. **The driver and the adapter live in `qtes-dftb`**, e.g. `tools/dftbplus/`
   for the milestone-1 driver.
2. **Target DFTB+ version: recommended `25.1`** (the release tag), recorded in
   the docs and the build notes.
   - The checkout examined is `25.1` plus 17 development commits; its library
     API is the same version (0.4.0, only two comments changed).
   - A tag gives a fixed, citable version that anyone can rebuild, matching
     the published manual and the conda-forge package.
   - The driver finds DFTB+ through pkg-config, so testing against a newer
     checkout only means pointing it to another install.
3. **SCC charges:** option A (§8.2) for milestone 1, to confirm the basic
   interface. No charges are written to or read from files. Setting the SCC
   start charges through the API (needed for options B and C) will be
   discussed with the DFTB+ developers.
4. **The charge-setting API call** will be developed by us on a fork branch
   and offered upstream as a pull request, after an issue agreeing on the
   design (§8.3).

**Still open:**

1. **Repeat counts in `.spl` files** (§6): does the DFTB+ SKF reader accept
   `19*0.0`? To be checked in milestone 1, step 3.
2. **Static or shared library:** static (the default) is simplest for a
   single executable; shared is only needed if Python or Julia should load
   it too.
3. **Charge option for production** (A, B or C), after the step 4
   measurements and the discussion with the DFTB+ developers.
