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
