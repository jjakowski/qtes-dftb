# agent_state.md — current state (overwrite on each update)

**Last updated:** 2026-11-17 (milestone 1 complete)
Check the latest commit with `git log --oneline -5`.

## Repository

- Path: `/home/j2c/Development/QTES/qtes-dftb`, branch `main`, remote
  `origin` (GitHub). Working tree should be clean after each commit.
- Environment: conda env `rmg-dev` (gfortran 15.2, OpenMPI, OpenBLAS,
  CMake 4.1.2, Python 3.12). Activate with
  `source ~/miniconda3/etc/profile.d/conda.sh && conda activate rmg-dev`.

## Status of the code

- Builds with gfortran + OpenBLAS (`(cd dftb && make libdftb) && make`).
- All examples have reference outputs; the regression test passes
  (byte-identical, except timing files).
- Bugs F1–F10 fixed; **open issues: #10** (sampling: a rejected point
  replaces the last accepted one; fixing it changes results, needs a
  decision) and **#11** (old DFTB engine crashed once with multithreaded
  OpenBLAS; always set `OPENBLAS_NUM_THREADS=1`). See
  `docs/KNOWN_ISSUES.md`.
- Docs: `README.md`, `docs/MANUAL.md` and `docs/manual/qtes-dftb-manual.tex`
  (draft; theory, parameter guidance, tutorials and developer notes still
  to be written with the user), `docs/DFTBPLUS_AUDIT.md`.

## Current project: use DFTB+ instead of the built-in DFTB code

Full analysis and decisions: `docs/DFTBPLUS_AUDIT.md` (§11 decisions).

### Decisions so far

1. Driver and adapter live in this repository (`tools/dftbplus/` for the
   milestone-1 driver).
2. Target DFTB+ **25.1** (release tag); found through pkg-config.
3. SCC charges: **option A** (DFTB+ default carry-over) for milestone 1; no
   charge files written or read.
4. An API call to set SCC start charges will be developed by us on a branch
   of the fork `jjakowski/dftbplus` and offered upstream (issue first).
5. The old engine stays selectable at first, for comparison.

### Roadmap

| # | Milestone | Status |
|---|---|---|
| 1 | Standalone driver: coordinates + settings in, energy + forces out from DFTB+ 25.1 | **done** (steps 1–4; step 5 optional, not done) |
| 2 | Adapter inside `qtes-dftb`: `derivs` calls DFTB+ via an input switch; settings in `IN_dftbplus.hsd` | not started |
| 3 | Dynamics with DFTB+ forces: examples, energy conservation, MPI behaviour, new references | not started |
| 4 | SCC start-charge API call in DFTB+ (issue → fork branch → PR) | not started (can start any time) |
| 5 | Charge handling in `qtes-dftb` (option B or C), compared with A | not started |
| 6 | Production and cleanup: performance, docs, keep or drop old engine | not started |

### Milestone 1 steps

1. **[DONE 2026-11-17] Build DFTB+ 25.1 library.** Exact commands, paths
   and results in `tools/dftbplus/BUILD.md`. Installed in
   `/home/j2c/Development/DFTB/dftbplus-25.1/install-qtes/`; `bin/dftb+`
   passed a CH4 SCC check. Original plan: export the tag
   with `git archive 25.1` from `/home/j2c/Development/DFTB/dftbplus` into
   e.g. `/home/j2c/Development/DFTB/dftbplus-25.1/`; build in
   `.../build-qtes/`, install to `.../install-qtes/`; `rmg-dev` gfortran,
   `WITH_API=ON WITH_MPI=OFF WITH_OMP=OFF`, OpenBLAS, static. Check that
   `libdftbplus.a`, module files and `dftbplus.pc` are installed and that
   the `dftb+` binary runs a DFTB+ test input. Record commands in
   `tools/dftbplus/BUILD.md`.
2. **[DONE 2026-11-17] Write the driver** `tools/dftbplus/dftbp_driver.f90`
   (+ `Makefile` via pkg-config, `IN_dftbplus.hsd`, `README.md`); builds and
   links with no warnings. Extras for step 3: writes the full DFTB+ input to
   `dftbp_input.hsd` (for check a) and has a finite-difference mode
   `dftbp_driver fd [h]` (for check b). Original plan: read geometry,
   `TDftbPlus_init` → `getInputFromFile` → add `Geometry` → `setupCalculator`
   → `setGeometry` → `getEnergy` / `getGradients`.
3. **[DONE 2026-11-17] First calculation and checks**: all four passed;
   results table in `tools/dftbplus/README.md`, runs in
   `tests/runs/dftbp-step3_20261117-181054/`. (a) driver = `dftb+` program
   to printing precision; (b) gradient vs FD 3e-9 Ha/bohr at h = 1e-4;
   (c) vs old engine (`tests/old_engine/`): E 3.4e-7 Ha, gradient 1.2e-5
   Ha/bohr; (d) `.spl` files read as they are. Rerun all checks with
   `tests/verify_dftbplus.sh` (PASS/FAIL). Original plan: first
   calculation on the `h-c4h9-14` geometry with the `.spl` files
   (6000 K, SCC, tolerance 1e-6), in `tests/runs/`: (a) same E and gradient
   as the `dftb+` binary on the same input; (b) gradient vs finite
   differences; (c) same order as the old engine; (d) DFTB+ accepts the
   `.spl` files (repeat counts like `19*0.0`).
4. **[DONE 2026-11-17] Repeated calls**: 300 `h-flake-111` geometries
   (111 atoms) in `tests/runs/dftbp-step4b_20261117-184152/`; results in
   `tools/dftbplus/README.md`. Option A vs fresh start: dE ≤ 1e-12 Ha,
   d grad ≤ 1.5e-7 Ha/bohr; same for call order. SCC iterations 14.3
   (fresh) / 11.9 (one instance, all QTs) / 8.0 (each QT in turn); time per
   call 0.25 / 0.20 / 0.14 s vs 0.42 s old engine. DFTB+ − old: constant
   3.16e-5 Ha offset, relative energies within 4.4e-7 Ha, gradients within
   2.2e-5 Ha/bohr. No files written. Original plan: several geometries on
   one instance (option A); SCC
   iterations and time per call vs old engine; order dependence within the
   SCC tolerance; no files written by DFTB+.
5. **Optional MPI version**: one instance per rank; no file clashes; what
   happens when DFTB+ stops on an error in one rank.
6. **Report and commit**: results in the audit report or
   `tools/dftbplus/README.md`; update README, `AGENTS.md` if rules change.

**Done when** step 3 (a) and (b) pass and step 4 shows stable repeated calls.

## Next action

Milestone 1 is complete (step 5, the optional MPI version of the driver,
was not done; its questions move to milestone 2/3). Next, after the user's
decision: **milestone 2**, the adapter inside `qtes-dftb` (`derivs` calls
DFTB+ through an input switch, settings in `IN_dftbplus.hsd`), and/or
milestone 4 (SCC start-charge API call in DFTB+, issue first).

Notes for milestone 2 (full list: `docs/DFTBPLUS_AUDIT.md` §12, lessons learned):
- Step-4 timing favours one DFTB+ instance per rank that follows its own
  QTs: with the contiguous QT blocks in `qtes-dftb`, each rank computes
  the same QTs every substep, i.e. close to the `byqt` case only if a rank
  holds one QT; with several QTs per rank it is like the `file` case.
- Always run with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` (#11).
- Driver checks: `tests/verify_dftbplus.sh`; series comparison:
  `tests/compare_dftbplus_series.py RUNDIR`.

## Open questions (for the user)

- Static or shared DFTB+ library? (static unless Python/Julia need it)
- Charge option for production: step 4 shows option A is accurate enough
  (effects ≤ 1.5e-7 Ha/bohr); B would mainly save SCC iterations. Decide
  together with the DFTB+ API discussion.
- Manual: meaning of the bonding probability `P`; is `ipot = 2` still used;
  guidance on QT count, `a0`, `dt`, `qmn`/`qmx`; theory chapter input.
- Issue #10 (sampling): fix, and accept changed results?

## Paths outside the repository (do not modify)

- DFTB+ checkout (upstream, clean): `/home/j2c/Development/DFTB/dftbplus`
  (`25.1-17-g7e1a071d`).
- Earlier DFTB+ work (forces-print feature, fork workflow notes):
  `/home/j2c/Development/DFTB/dftbplus_forces_print`,
  `/home/j2c/Development/DFTB/readme-short.md`.
- DFTB+ 25.1 source/build/install for this project (created in step 1):
  `/home/j2c/Development/DFTB/dftbplus-25.1/` (`build-qtes/`,
  `install-qtes/`, build logs `*-qtes.log`).
- Parameter set used for checks: `/home/j2c/Development/DFTB/SlaKo/3ob-3-1/`.
