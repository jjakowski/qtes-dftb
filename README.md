# qtes-dftb

`qtes-dftb` is a parallel quantum-trajectory dynamics code (`bohmd`, version
17/08/2016) with forces from DFTB. Some atoms are treated as quantum
particles. Their wavepacket is sampled by an ensemble of quantum trajectories
(QTs), and the remaining atoms move classically. In the Ehrenfest mode, the
classical atoms feel the force averaged over the QTs. All the included
examples simulate one H atom scattering off a hydrocarbon (C/H) target.

The code is Fortran 77/90 and uses MPI over trajectories, OpenMP, and
BLAS/LAPACK.

## Repository layout

```
qtes.f              main program: sampling, propagation, quantum potential, output
read_input_mpi.f    reads INdyn / INxyz on rank 0, sets up the restart from traj.sav
set_traj.f          reloads trajectories from traj.sav on restart
dumpef.f            debug dump of energies and forces (fort.1000+rank)
bohm.inc            array limits: ndx=700 dims (~233 atoms), npx=5000 QTs, nqx=200 quantum dims
Makefile            builds qtes.x
dftb/               DFTB force library (libdftb.a) with its own GNUmakefile / GNUmake.inc
examples/           input sets and reference results, see examples/README.md
tests/              run_example.sh; runs go to tests/runs/ (git-ignored)
docs/               DYNAMICS.md, KNOWN_ISSUES.md; original notes: file formats (sg-notes-files), code notes, sbatch script
```

## Requirements

- An MPI Fortran compiler wrapper (`mpif90`). The build is tested with
  OpenMPI and gfortran 15.
- OpenBLAS, which provides both BLAS and LAPACK.

On the development machine, all of these come from the `rmg-dev` conda
environment. Activate it before you build or run, because `qtes.x` loads
OpenBLAS and the MPI libraries from that environment at run time.

## Building

```bash
(cd dftb && make libdftb)   # builds dftb/libdftb.a; the top-level make does not build it
make                        # builds qtes.x
```

The build uses `-fallow-argument-mismatch -std=legacy`. gfortran 10 and
later reject the legacy argument type and rank mismatches in this code unless
these flags are set. The compiler prints many warnings, which is expected.
The original Intel/MKL settings (`-mkl=sequential`, `-openmp`) are kept as
comments in `Makefile` and `dftb/GNUmake.inc`.

To rebuild from scratch:

```bash
make allclean && (cd dftb && make clean && rm -f *.o libdftb.a)
```

Note that `make allclean` also removes `fort.*` and `*.log` from the top-level
directory.

## Running

The program takes no command-line arguments. It reads fixed file names
(`INdyn`, `INxyz`, `IN_dftb`, and the `*.spl` Slater-Koster files) from the
current directory and writes all output there as well. To keep the reference
data safe, always run through the test script:

```bash
tests/run_example.sh <example> [nprocs]     # e.g. tests/run_example.sh h-c4h9-14-ref 4
```

The script creates `tests/runs/<example>_<timestamp>/`, copies in the
example's inputs, the SK files, any restart checkpoint, and `qtes.x`, and then
runs `mpirun` there. It also sets `OMP_NUM_THREADS=1` (unless you set it
already), sets `OPENBLAS_NUM_THREADS=1`, and raises the stack limit with
`ulimit -s unlimited`.

To set up a run by hand, copy the same files into an empty directory and run
`mpirun -np N ./qtes.x > stdout.log`. Use the `mpirun` from the same
environment that provided `mpif90`.

### Restarting a run

To restart a run, set the restart flag in `INdyn` (line 9, `<nQT> 1`). The
code then reads the time and all trajectories from `traj.sav` and runs the
remaining `kmax - t/dt` large steps. It appends results to `traj`, `aver`,
`ener`, `cor`, `qt.dat`, `qtx`...`clz`, `pfit`, and `wft`, so those files
must already exist in the run directory. See the `h-c4h9-14-restart`
example.

## Examples and regression tests

| Example | Purpose |
|---|---|
| `h-c4h9-14-ref` | main regression test, fresh run from t = 0 to 100, full reference output |
| `h-c4h9-14-restart` | restart test, t = 80 to 100 from a checkpoint; matches the fresh run |
| `h-c4h9-14` | quick smoke test (about 1 s) |
| `h-flake-53` | 53-atom system, full reference output (run on 48 ranks) |
| `h-flake-53-print` | same system, demonstrates the `nqtprt`/`nclprt` print options (prints all QTs and all classical atoms) |
| `h-flake-111` | chemisorbed H on a 110-atom flake, quantum potential + sampling (`ipot=1`, `isr=1`), full reference output (run on 12 ranks) |

`examples/README.md` gives the details for each example. A quick check after a
code change:

```bash
tests/run_example.sh h-c4h9-14-ref 4
diff -q tests/runs/h-c4h9-14-ref_*/ener examples/h-c4h9-14-ref/reference/ener
```

The current reference results come from the gfortran + OpenBLAS build. They
agree with the original Intel/MKL results (commit `2ed84f9`) to round-off.
The timing files `trajt` and `dftbt` always differ from run to run.

## Output files

The main outputs are `ener` (energies), `cor` (correlation function), `aver`
(norm, bonding probability, average positions), `traj` (representative QT),
`qt.dat` (all QTs at every large step), `traj.sav` (restart checkpoint),
`wf0` and `wft` (initial and final QT positions), and `qtx`...`clz`.
`docs/sg-notes-files` gives the column-by-column formats.

- `qtx/qty/qtz` give the quantum atom's x/y/z coordinate for up to
  `nqtprt` QTs.
- `clx/cly/clz` depend on the mode. In Ehrenfest mode (`ipf = 0`), they
  give up to `nclprt` classical atoms from one QT, because the classical
  atoms are the same in every QT. Otherwise, they give the first classical
  atom for up to `nqtprt` QTs.

Both counts are optional values on line 9 of `INdyn`, after the number of
QTs and the restart flag:

```
48 0              \  defaults: nqtprt = 20, nclprt = 20
48 0 48 0         \  48 QTs in qtx..qtz, all classical atoms in clx..clz
```

A value of 0 or less means "all". The values in use are printed in
`stdout.log` as `PRINT OPTIONS`. For a full usage guide, see "Print
options" in `examples/README.md`, which includes the column mapping, the
restart note, and the `h-flake-53-print` example to try.

## Dynamics modes

[`docs/DYNAMICS.md`](docs/DYNAMICS.md) describes the run-mode switches and
what the code does for each one:
- `ipf`: Ehrenfest mode vs independent classical copies
- `ipot`: the quantum potential
- `isr` / `pow`: initial sampling
- `qmn` / `qmx`: the z-window that excludes trajectories
- the restart flag

## Known issues

Open bugs and behaviour that is easy to miss are listed in
[`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md).

## Known messages

- gfortran reports at exit any floating-point exceptions that were raised
  during the run.
  - The former `IEEE_DIVIDE_BY_ZERO` / `IEEE_INVALID_FLAG` messages are
    fixed (`docs/KNOWN_ISSUES.md` F4).
  - `IEEE_DENORMAL`, seen in the restart example, is harmless.
  - To locate a new exception, build with
    `-g -fbacktrace -ffpe-trap=zero,invalid`.
- `MPI_FINALIZE ERR 0` is the return code of `MPI_Finalize`, printed once
  per rank. The value 0 means success.
