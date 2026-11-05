# Examples

Every example is a scattering simulation of one quantum H atom hitting a
hydrocarbon target. All other atoms are treated classically, and DFTB supplies
the forces. Run the examples with `tests/run_example.sh` (see below). Do not
run anything inside `examples/`: the code writes output files to the directory
it runs in, which would overwrite the reference data.

## Layout

```
examples/
  sk/                 shared Slater-Koster tables (*.spl)
  <example>/
    INdyn             dynamics parameters + initial wavefunction (one line per atom)
    INxyz             model-potential parameters + starting geometry (Angstrom)
    IN_dftb           DFTB setup: ntype, lmax, which SK files to use
    reference/        results from earlier runs (only where available)
```

The code expects the input files to have exactly these names. `IN_dftb`
selects `cc.spl`, `ch-g2ms.spl`, `hc-g2ms.spl` and `hh.spl`. The `*-b3.spl`
files and the plain `ch.spl`/`hc.spl` are alternative C-H parameterizations.
To use one of them, edit `IN_dftb`.

`docs/sg-notes-files` describes every input and output file column by column.

## Examples

| Example | System | Steps (kmax x kout, dt) | QTs | Reference |
|---|---|---|---|---|
| `h-c4h9-14-ref` | H + C4H9 cluster, 14 atoms | 5 x 10, dt 20 (t = 0 to 100) | 48 | full output set |
| `h-c4h9-14`     | same system | 2 x 4, dt 10 (short smoke test) | 48 | none |
| `h-flake-53`    | H + 52-atom C/H flake | 10 x 3, dt 10 | 48 | screen output only |
| `h-flake-111`   | H + 110-atom C/H flake, `ipot=1` (linear p-fit QP) | 25 x 2, dt 3 | 12 | none |

- **h-c4h9-14-ref**: this is the main regression test. `reference/` holds
  the full output set from the original run (`ener`, `cor`, `aver`, `traj`,
  `qt.dat`, `traj.sav`, `wft`, ...). The original `INdyn` had the restart flag
  set to 1. It is now set to 0, so a run starts from t = 0 and reproduces the
  reference. The original data was produced in two stages: a run from t = 0
  to 80, then a restart from 80 to 100. As a result, `reference/wf0` holds
  the t = 80 restart state, not t = 0. `reference/wft` holds two 48-row
  blocks, t = 80 then t = 100, and a fresh run matches only the second
  block. Also, `qtx/qty/qtz` and `clx/cly/clz` contain only the time
  column, because of a bug in `aver` in `qtes.f`. The current gfortran +
  OpenBLAS build reproduces `ener`, `aver`, `traj`, `cor` and `qt.dat` to
  round-off (max relative difference about 1e-10).
- **h-c4h9-14**: uses the same geometry as above with a much shorter run
  (about 1 s on 4 ranks). Use it for a quick build check.
- **h-flake-53**: `reference/stdout.log` is the saved screen output of an
  earlier run (formerly `junk`). It contains no output files.
- **h-flake-111**: the largest system. This is the only example that uses a
  nonzero quantum potential (`ipot=1`).

## Running

Build first (from the repository root):

```bash
(cd dftb && make libdftb) && make
```

Then run, for example:

```bash
tests/run_example.sh h-c4h9-14-ref 4     # <example> [MPI ranks, default 4]
```

Each run gets a new directory, `tests/runs/<example>_<timestamp>/`, which git
ignores. To compare a run with the reference:

```bash
diff tests/runs/h-c4h9-14-ref_*/ener examples/h-c4h9-14-ref/reference/ener
```
