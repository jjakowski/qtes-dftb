# CLAUDE.md

Instructions for Claude Code in this repository. User-facing documentation
is in `README.md` and the manual (`docs/MANUAL.md`, `docs/manual/`).

## Working rules

- **Ask before acting.** Explain the plan and wait for approval before
  changing code, fixing a bug, committing, or pushing. Approval of one step
  does not cover the next.
- **No attribution.** Never add `Co-Authored-By: Claude` or any other
  Claude/Anthropic attribution to commits or pull requests.
- **Bugs go in `docs/KNOWN_ISSUES.md`.** Record each bug found under
  **Open**. Fix it only with approval, then move it to **Fixed** (F1, F2,
  ...) with the commit and how it was verified.
- **Bug fixes must not change the physics** unless that is the intent.
  Show this with a regression test (below) before asking to commit.
- **Keep the docs correct.** When code layout, input parameters or features
  change, update in the same commit everything that describes them:
  `README.md`, `docs/MANUAL.md`, `docs/manual/qtes-dftb-manual.tex`, and
  where relevant `examples/README.md` and `docs/KNOWN_ISSUES.md`. The README
  may be shorter than the manual, but must not be wrong.

## Build

All tools come from the conda environment `rmg-dev` (gfortran 15, OpenMPI,
OpenBLAS). Activate it before building and before running.

```bash
(cd dftb && make libdftb)   # top-level make does not build this
make                        # -> qtes.x
```

Flags: `-O2 -fopenmp -fallow-argument-mismatch -std=legacy`, linked with
`-lopenblas`. Many compiler warnings are expected. For floating-point bugs,
build a separate copy with `-O0 -g -fbacktrace -ffpe-trap=zero,invalid`.

Manual PDF: `cd docs/manual && make` (the PDF and LaTeX build files are
git-ignored).

## Running and testing

- **All test and debug runs go in `tests/`.** Never run `qtes.x` inside
  `examples/`: it writes output into the current directory and would
  overwrite reference data.
- `tests/run_example.sh <example> [ranks]` creates
  `tests/runs/<example>_<timestamp>/` (git-ignored) and runs there.
- **Rank count:** every MPI rank must get at least one QT, or the run hangs
  (open issue #8). Use a divisor of the QT count.
- **Regression test:** run every example and compare with
  `examples/<example>/reference/`. All files must be byte-identical except
  the timing files `trajt` and `dftbt`. The loop is in the manual, chapter
  "Verification and regression tests". Reference rank counts:
  `h-c4h9-14-ref 4`, `h-c4h9-14-restart 4`, `h-flake-53 48`,
  `h-flake-53-print 48`, `h-flake-111 12`.
- For changes that are meant to alter results, check energy conservation
  (`E_total` in `ener`) instead, and replace references only with approval.

## Code conventions

- Fixed-form Fortran 77: code past column 72 is silently ignored. Keep new
  lines short, and match the style of the surrounding code.
- Implicit typing is in effect: a misspelled or undeclared variable is not
  an error. Check new variable names carefully.
- MPI uses replicated data: every rank holds all QTs, and DFTB work is split
  into contiguous blocks. Any decision that changes the state (weights,
  exclusion, positions) must be made identically on every rank (see F7).
- Shared settings live in common blocks; when adding one, declare it in
  every routine that uses it and broadcast inputs from rank 0.

## Shell note

`grep` is aliased to `ugrep`, which skips files it considers binary (some
output files). Use `command grep -a`.
