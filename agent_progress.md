# agent_progress.md — history (append-only)

Each entry: date, commit(s), what changed and why, tests and results, next
steps. Newest entries at the bottom. Details of each commit are in
`git log`; bug details are in `docs/KNOWN_ISSUES.md`.

## 2026-10-30 — Build without MKL

- Commit `43cbed5`. Switched from Intel/MKL to gfortran + OpenBLAS + OpenMPI
  (`rmg-dev`). Needed `-fallow-argument-mismatch -std=legacy` for legacy
  argument mismatches.
- Test: first runs of the original example inputs.

## 2026-11-05 — Examples, tests, bug fixes, print options

- Commits `ab98622` … `b851640`.
- Reorganized into `examples/` (inputs + `reference/`) and `tests/`
  (`run_example.sh`, runs in `tests/runs/`). Added README files,
  `docs/KNOWN_ISSUES.md`, `docs/DYNAMICS.md` (later folded into the manual).
- Examples: `h-c4h9-14-ref`, `h-c4h9-14-restart` (restart from t = 80),
  `h-c4h9-14` (smoke test), `h-flake-53`, `h-flake-53-print`, `h-flake-111`.
- Fixed F1 (`aver` printed no trajectories), F2 (`h-flake-111` inputs), F3
  (undefined weight index in `trot`), F4 (divide-by-zero in `init_traj`), F5
  (forces outside the `qmn`/`qmx` window), F6 (fit failure: clean stop on all
  ranks), F7 (exclusion below `qmn` on every rank), F8 (quantum-atom force
  applied twice with `ipot = 1`).
- Added optional print counts `nqtprt`, `nclprt` on `INdyn` line 9.
- Tests: regression byte-identical where physics unchanged; energy
  conservation for F8; per-rank checksums for F7.

## 2026-11-12 … 2026-11-17 — Manual, CLAUDE.md

- Commits `c57604b`, `e61640c`, `5ce6a81`, `8d20d2d`.
- User manual draft in Markdown (`docs/MANUAL.md`) and LaTeX
  (`docs/manual/`, PDF via `make`); run modes chapter and Known limitations
  appendix; `CLAUDE.md` with working rules.
- Citations: Garashchuk et al., JCTC 9, 5221 (2013); Wang et al., JPCC 118,
  16175 (2014).

## 2026-11-17 — Issue #8 (F9): empty MPI rank

- Commit `5f2c773`. New `rankchk`, `qtsplit`, `prange` in `qtes.f`: stop
  with `RANK ERROR` and suggestions (rank counts, QT counts) when a rank
  would get no QT; one-line `QT SPLIT` in every `stdout.log`.
- Tests: former hangs now stop cleanly; brute-force check of the
  suggestions (908,658 cases); regression byte-identical.

## 2026-11-17 — Issue #9 (F10): restart file checks; issue #10 recorded

- Commit `54a3315`. `traj.sav` now starts with `# np nd dt kout`; shared
  reader `rdsav` (in `set_traj.f`) checks it against the inputs, requires
  exactly `np` complete QTs, stops all ranks cleanly. Old files still read.
- Tests: all restart error cases; regression byte-identical except the
  `traj.sav` header (references replaced).
- Recorded open issue #10 (sampling).

## 2026-11-17 — DFTB+ integration audit

- Commits `d0c7c12`, `496e099`. `docs/DFTBPLUS_AUDIT.md`: DFTB+ 25.1 API fits
  (`setGeometry`, `getEnergy` = Mermin, `getGradients`, bohr/hartree);
  decisions recorded (driver in `tools/dftbplus/`, target 25.1, charge option
  A, own API call for SCC start charges via fork + PR); Fortran 77/90 notes.
- No code changed, nothing built.

## 2026-11-17 — Agent handoff files

- Added `AGENTS.md` (rules for any agent; content moved from `CLAUDE.md`),
  `CLAUDE.md` now points to it, `agent_state.md`, `agent_progress.md`,
  `docs/RESTART_PROMPT.md`.
- Next: milestone 1, step 1 (DFTB+ 25.1 library build), after approval.

## 2026-11-17 — Milestone 1, step 1: DFTB+ 25.1 library built

- Exported tag `25.1` (commit `8802aff6`) with `git archive` into
  `/home/j2c/Development/DFTB/dftbplus-25.1/`; the checkout stays clean.
- Configured with `rmg-dev` gfortran 15.2, `WITH_API=ON WITH_MPI=OFF
  WITH_OMP=OFF`, static, OpenBLAS; built in about 1 min (`-j 32`) with no
  errors; installed to `install-qtes/` (`libdftbplus.a`, 299 module files,
  `dftbplus.h`, `dftbplus.pc` version 25.1, `bin/dftb+`).
- Check: CH4 SCC single point with `bin/dftb+` and 3ob-3-1 C/H (Fermi
  6000 K, tolerance 1e-6) in `tests/runs/dftbp-step1-ch4_*`: converged in 4
  iterations, E_Mermin = -3.2291223287 Ha, forces tetrahedrally symmetric.
  First attempt failed on my input keyword (`Analysis/CalculateForces`);
  fixed to `PrintForces`.
- Added `tools/dftbplus/BUILD.md` with exact commands and results.
- Next: step 2, the driver.

## 2026-11-17 — Milestone 1, step 2: driver written, builds and links

- Committed step 1 as `c1d2e76`.
- New in `tools/dftbplus/`: `dftbp_driver.f90` (free-form; reads `INxyz` and
  `IN_dftbplus.hsd`, adds the Geometry block in code, `setupCalculator`,
  then `setGeometry`/`getEnergy`/`getGradients`; second geometry with atom 1
  moved; optional finite-difference mode; dumps the full input to
  `dftbp_input.hsd`), `IN_dftbplus.hsd` (old hard-coded settings as input:
  SCC, tolerance 1e-6, charge 0, Fermi 6000 K, max. angular momenta, `.spl`
  files; DFTB+ output files switched off), `Makefile` (pkg-config),
  `README.md`.
- Build: `make` in `tools/dftbplus` with rmg-dev gfortran: no errors or
  warnings. First attempt failed because `export PKG_CONFIG_PATH` in the
  Makefile does not reach `$(shell ...)`; fixed by passing it in the command.
- `.gitignore`: `tools/dftbplus/dftbp_driver`, `*.mod`.
- Not yet run (that is step 3).
- Next: step 3, the checks.

## 2026-11-17 — Milestone 1, step 3: driver checks, all passed

- Committed step 2 as `72f2a96`.
- Runs in `tests/runs/dftbp-step3_20261117-181054/` (h-c4h9-14 `INxyz`,
  `examples/sk/*.spl`, `IN_dftbplus.hsd`).
- Problems found and fixed on the way:
  1. Driver read all of `INxyz`; that file has a second geometry copy after
     a blank line and a comment. Now stops at the first blank line after
     the atoms (as `qtes.x`) and checks `nat` from `INdyn` if present.
  2. `getGradients` stopped: "Forces not available". In DFTB+ 25.1 forces
     without a Driver are switched on by `Analysis { PrintForces = Yes }`
     (`initprogram.F90`, `tForces = ctrl%tForces .or. tPrintForces`);
     `Options { CalculateForces }` (used in DFTB+'s API tests with
     ParserVersion 3) is rejected with ParserVersion 14. Added to
     `IN_dftbplus.hsd`. The step-1 note claiming CalculateForces is not a
     keyword was inaccurate (it lives in old syntax only); corrected.
- Results: (d) `.spl` files read as they are, SCC in 8 iterations;
  (a) driver vs `bin/dftb+` on `dftbp_input.hsd`: E within 3e-11 Ha,
  gradient within 5e-13 Ha/bohr; (b) FD: 2.3e-7 / 3.1e-9 / 1.5e-9 Ha/bohr
  for h = 1e-3 / 1e-4 / 1e-5; (c) new `tests/old_engine/` (calls `do_dftb`
  as `derivs` does): E −10.708373138 vs −10.708373482 Ha, energy change
  1.3030151e-3 vs 1.3030623e-3 Ha, gradient max diff 1.2e-5 Ha/bohr.
- Docs: results in `tools/dftbplus/README.md`; `tests/old_engine/README.md`;
  `BUILD.md` note on forces; audit §11 question on `.spl` answered;
  README layout.
- Next: step 4 (repeated calls, timing, order dependence).
- Added `tests/verify_dftbplus.sh`: builds both programs, runs checks a–d
  in `tests/runs/verify-dftbplus_<timestamp>/` copying only the files each
  check needs, prints PASS/FAIL with fixed tolerances. First version failed
  because `set -u` breaks `conda activate`; moved after activation. Result:
  all four PASS, same numbers as above.
- Confirmed DFTB+ reads both `hc-g2ms.spl` and `ch-g2ms.spl` (the log lists
  only the first file of each pair; a run without `ch-g2ms.spl` stops with
  "SK file not found").

## 2026-11-17 — Milestone 1, step 4: repeated calls; milestone 1 complete

- Committed step 3 as `063c31b`.
- Driver: new `series` mode (read a `qt.dat`, orders `file`/`byqt`/
  `shuffle`, optional `fresh` instance per geometry; writes
  `series_result.dat` with time per call and `series_grad.dat`); instance
  setup moved into `setupInstance`. `tests/old_engine`: `series` mode.
  New `tests/compare_dftbplus_series.py` (energies, gradients, SCC
  iterations from the DFTB+ log, timings).
- Problems on the way:
  1. First parallel runs forgot `OPENBLAS_NUM_THREADS=1`: jobs up to 100×
     slower, and the old engine crashed with SIGSEGV (`tests/runs/
     dftbp-step4_20261117-182718/old/`). With the variable set it runs
     fine. Recorded as open issue #11 (cause not established); manuals and
     `verify_dftbplus.sh` updated.
  2. My first rerun command backgrounded the whole setup chain with `&`,
     so only one job got the folder variable; fixed with a script.
- Results (300 `h-flake-111` geometries, `tests/runs/dftbp-step4b_
  20261117-184152/`): option A vs fresh ≤ 1.1e-12 Ha / 1.5e-7 Ha/bohr;
  call order ≤ 1.0e-12 Ha / 1.4e-7 Ha/bohr; SCC iterations 14.3 / 11.9 /
  11.7 / 8.0 (fresh / file / shuffle / byqt); time per call 0.251 / 0.201 /
  0.202 / 0.141 s, old engine 0.420 s; DFTB+ − old: constant 3.16e-5 Ha
  (spread 1.4e-7), relative energies 4.4e-7 Ha, gradients 2.2e-5 Ha/bohr.
  No DFTB+ files written.
- `tests/verify_dftbplus.sh` still passes after the driver changes.
- Docs: step-4 section in `tools/dftbplus/README.md`, old-engine README,
  KNOWN_ISSUES #11 + manuals (threads note, appendix D), audit §11, README.
- Next: milestone 2 or 4, user to decide.
- Added §12 "Lessons learned from milestone 1" to `docs/DFTBPLUS_AUDIT.md`.
