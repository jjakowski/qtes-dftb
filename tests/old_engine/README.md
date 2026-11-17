# tests/old_engine: the old DFTB engine on one geometry

Test program for milestone 1, step 3, check (c): computes the energy and
gradient with the old built-in DFTB code (`dftb/libdftb.a`), set up and
called exactly as `qtes.x` does it (`dftb_init_mpi` as in the `qtes.f` main
program, `do_dftb` with the settings of `derivs`: 6000 K, SCC tolerance
1e-6, neutral). `qtes.f` is not changed.

It prints the same two geometries as `tools/dftbplus/dftbp_driver` (the
`INxyz` geometry, then atom 1 moved by +0.1 bohr in z) and writes
`old_result.dat` in the same format as `dftbp_result.dat`.

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate rmg-dev
cd tests/old_engine && make            # needs dftb/libdftb.a
R=tests/runs/old_engine_$(date +%Y%m%d-%H%M%S); mkdir -p $R
cp examples/h-c4h9-14/{INxyz,INdyn,IN_dftb} examples/sk/*.spl $R/
cd $R && ulimit -s unlimited && mpirun -np 1 ../../old_engine/old_engine
```

`old_engine series qt.dat` instead computes every geometry of a
`qtes-dftb` `qt.dat` (milestone 1, step 4) and writes
`old_series_result.dat` (index, energy, time per call) and
`old_series_grad.dat`. Set `OPENBLAS_NUM_THREADS=1` (see
`docs/KNOWN_ISSUES.md` #11).

Inputs: `INdyn` (line 1: `nat`), `INxyz` (first `nat` atoms, Angstrom,
converted with 0.52918 Å/bohr as in `read_input_mpi.f`), `IN_dftb` and the
`.spl` files it names. Results of the comparison: `tools/dftbplus/README.md`.
