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
    restart/          restart examples only: checkpoint traj.sav + output files to append to
    reference/        results from earlier runs (only where available)
```

Line 9 of `INdyn` can also hold two optional print counts, `nqtprt` and
`nclprt`. See "Print options" below. Only `h-flake-53-print` sets them;
the other examples use the defaults.

The code expects the input files to have exactly these names. `IN_dftb`
selects `cc.spl`, `ch-g2ms.spl`, `hc-g2ms.spl` and `hh.spl`. The `*-b3.spl`
files and the plain `ch.spl`/`hc.spl` are alternative C-H parameterizations.
To use one of them, edit `IN_dftb`.

`docs/sg-notes-files` describes every input and output file column by column.

## Examples

| Example | System | Steps (kmax x kout, dt) | QTs | Reference |
|---|---|---|---|---|
| `h-c4h9-14-ref` | H + C4H9 cluster, 14 atoms | 5 x 10, dt 20 (t = 0 to 100) | 48 | full output set |
| `h-c4h9-14-restart` | same system, restart from t = 80 | 5 x 10, dt 20 (runs t = 80 to 100) | 48 | full output set |
| `h-c4h9-14`     | same system | 2 x 4, dt 10 (short smoke test) | 48 | none |
| `h-flake-53`    | H + 52-atom C/H flake | 10 x 3, dt 10 (t = 0 to 100) | 48 | full output set + old screen log |
| `h-flake-53-print` | same as `h-flake-53`, with print options `nqtprt=48`, `nclprt=0` (all) | 10 x 3, dt 10 (t = 0 to 100) | 48 | full output set |
| `h-flake-111`   | chemisorbed H on a 110-atom C/H flake, `ipot=1`, `isr=1` | 25 x 2, dt 3 (t = 0 to 75) | 12 | full output set |

- **h-c4h9-14-ref**: this is the main regression test. `reference/` holds
  the full output set (`ener`, `cor`, `aver`, `traj`, `qt.dat`, `traj.sav`,
  `wf0`, `wft`, `qtx`...`clz`, plus `stdout.log`) from a fresh t = 0 to 100
  run on 4 MPI ranks. The run used the gfortran + OpenBLAS build after the
  `aver` fix, which restored the trajectory columns in `qtx/qty/qtz` and
  `clx/cly/clz`. These results agree with the original Intel/MKL results to
  round-off (max relative difference about 1e-10). The original results are
  kept in git history, in commit `2ed84f9`. They were produced in two
  stages, a run from t = 0 to 80 and a restart from 80 to 100, so their
  `wf0` and `wft` differ in form. The timing files `trajt` and `dftbt`
  differ from run to run.
- **h-c4h9-14-restart**: tests the restart path (`set_traj`, reading
  `traj.sav`). The inputs are the same as `h-c4h9-14-ref`, except that the
  restart flag in `INdyn` is set to 1 (`48 1`). `restart/` holds a
  checkpoint from a fresh run of t = 0 to 80: `traj.sav` at t = 80 and the
  stage-1 output files. A restart needs those output files, because the
  code opens them with `status='old'` and appends to them. The code reads
  t = 80 from `traj.sav` and runs the remaining `kmax - t/dt` = 1 large
  step, which is 10 substeps. After the run, `ener`, `aver`, `traj`, `cor`,
  `qt.dat`, `traj.sav` and `qtx`...`clz` cover t = 0 to 100. The rows that
  the restart job itself writes (t > 80) were checked against the single
  fresh 0 to 100 run in `h-c4h9-14-ref/reference/`. `ener` and `cor` are
  bit-identical, and the final state in `traj.sav` agrees to 1e-13 or
  better. Two files differ by design. `wf0` is the t = 80 restart state,
  because `set_traj` rewrites it. `wft` has a t = 80 block followed by a
  t = 100 block, because it is appended.
- **h-c4h9-14**: uses the same geometry as above with a much shorter run
  (about 1 s on 4 ranks). Use it for a quick build check.
- **h-flake-53**: `reference/` holds the full output set from a fresh run
  on 48 MPI ranks (about 4 s) with the gfortran + OpenBLAS build.
  `reference/stdout_old.log` (formerly `junk`) is the screen output of the
  original run, also on 48 ranks, from a slightly older code version that
  printed extra debug dumps. That log has no output files, so only the
  t = 0 quantities can be compared with it. Those quantities agree: the
  sampling, the norm and the kinetic energy are identical, and the DFTB
  potential energy agrees to 2e-12 relative. Over the run, the total
  energy in `ener` stays between 0.050685 and 0.050706.
- **h-flake-53-print**: demonstrates the print options. The inputs are
  identical to `h-flake-53`, except for line 9 of `INdyn`: `48 0 48 0`,
  which prints all 48 QTs in `qtx..qtz` and all 52 classical atoms in
  `clx..clz`. As a result, `qtx..qtz` have 49 columns instead of 21, and
  `clx..clz` have 53 instead of 21. All other output files are
  byte-identical to the `h-flake-53` reference, and so are the first 21
  columns of `qtx..clz`: the options change only what is printed. Run it
  on 48 ranks (about 4 s).
- **h-flake-111**: the largest system. Unlike the other examples, the H
  is not incoming. It starts chemisorbed on C2 (C-H about 1.12 A) at
  rest. This is the only example that uses the quantum potential
  (`ipot=1`, a linear quantum-force fit) and Gaussian sampling
  (`isr=1`), with 12 QTs spread around the H position. A run takes
  about 25 s on 12 ranks: the C-H bond vibrates, and the total energy
  is conserved to all 7 printed digits. The input was repaired on
  2026-11-05 (`docs/KNOWN_ISSUES.md` F2). `reference/` was created after
  the `ipot=1` double-force fix (F8), so use 12 MPI ranks to compare
  with it.

## Print options (`nqtprt`, `nclprt`)

Two optional values on line 9 of `INdyn`, after the number of QTs and the
restart flag, set how many columns the trajectory snapshot files get:

```
 np  restart  [nqtprt  [nclprt]]   \ comment
 48  0                             \ defaults: nqtprt = 20, nclprt = 20
 48  0        48       0           \ all 48 QTs, all classical atoms
```

| Value | Controls | Default | `<= 0` |
|---|---|---|---|
| `nqtprt` | QTs in `qtx/qty/qtz` (and in `clx/cly/clz` when `ipf != 0`) | 20 | all QTs |
| `nclprt` | classical atoms in `clx/cly/clz` (Ehrenfest mode, `ipf = 0`) | 20 | all classical atoms |

- **Order:** to set `nclprt`, you must also give `nqtprt`. Values larger
  than the real number are capped automatically.
- **Check:** `stdout.log` prints the values in use, for example
  `PRINT OPTIONS nqtprt, nclprt (<=0: all): 48 0`.
- **`qtx/qty/qtz` columns:** time, then x/y/z of the quantum atom in QT 1,
  QT 2, ...
- **`clx/cly/clz` columns** (`ipf = 0`): time, then x/y/z of classical
  atoms 1, 2, ..., all taken from the first QT with nonzero weight.
  Classical atom *m* is atom number *nq + m* in the `INdyn`/`INxyz` atom
  list. With one quantum atom, column 2 is atom 2, column 3 is atom 3,
  and so on.
- **`clx/cly/clz` with `ipf != 0`:** time, then the first classical atom
  in QT 1, QT 2, ... (`nqtprt` columns).
- **Restart:** a restart appends rows to `qtx..clz`. Use the same
  values in the original run and in the restart, or the rows will differ
  in length.
- **Not affected:** `wf0`/`wft` and `aver` keep their own limit of 40
  coordinates. `qt.dat` and `traj.sav` always contain everything.

Try it with the example:

```bash
tests/run_example.sh h-flake-53-print 48
R=$(ls -d tests/runs/h-flake-53-print_* | tail -1)
grep "PRINT OPTIONS" $R/stdout.log
awk '{print NF}' $R/qtx | sort -u      # 49 = time + 48 QTs
awk '{print NF}' $R/clx | sort -u      # 53 = time + 52 classical atoms
```

To use other values, copy an example to a work directory and edit line 9
of its `INdyn`. Do not edit the files in `examples/`: the references were
made with these inputs.

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
