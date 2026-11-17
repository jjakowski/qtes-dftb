# tools/dftbplus: calling DFTB+ from a qtes-dftb-style driver

Milestone 1 of the DFTB+ project (`docs/DFTBPLUS_AUDIT.md`,
`agent_state.md`): show that a small Fortran program can pass coordinates
and settings to DFTB+ and get the energy and forces back.

| File | Content |
|---|---|
| `BUILD.md` | how the DFTB+ 25.1 library was built and installed (outside this repository) |
| `dftbp_driver.f90` | the driver |
| `IN_dftbplus.hsd` | DFTB+ settings (template; copy into the run directory) |
| `Makefile` | builds `dftbp_driver` against the installed library via pkg-config |

## Build

```bash
source ~/miniconda3/etc/profile.d/conda.sh && conda activate rmg-dev
cd tools/dftbplus && make          # -> dftbp_driver (git-ignored)
```

Another DFTB+ install can be used with
`make DFTBPLUS_PREFIX=/path/to/install`.

## Run

Run in a directory under `tests/runs/` that contains `INxyz` (qtes-dftb
format, Angstrom), `IN_dftbplus.hsd`, and the Slater–Koster files it names
(e.g. from `examples/sk/`). Optionally `INdyn`: its line 1 (`nat`) is
checked against the atom count. As in `qtes.x`, the geometry in `INxyz`
ends at the first blank line after the atoms; anything after it is ignored.

```bash
R=tests/runs/dftbp_$(date +%Y%m%d-%H%M%S); mkdir -p $R
cp examples/h-c4h9-14/INxyz tools/dftbplus/IN_dftbplus.hsd examples/sk/*.spl $R/
cd $R && ulimit -s unlimited && ../../../tools/dftbplus/dftbp_driver
```

Modes:

```bash
dftbp_driver          # energy + gradient at the INxyz geometry, then with atom 1 moved +0.1 bohr in z
dftbp_driver fd       # also compares the gradient with finite differences (step 1e-4 bohr)
dftbp_driver fd 1e-3  # same, with step 1e-3 bohr
```

The geometry is taken from `INxyz` as is; the `INdyn` displacements used by
`qtes.x` are not applied.

## What the driver does

1. Reads `INxyz` and converts to bohr (0.52918 Å/bohr, as `qtes-dftb`).
2. `TDftbPlus_init` (DFTB+ screen output goes to `dftbp.log`).
3. Reads `IN_dftbplus.hsd` (`getInputFromFile`) and adds the `Geometry`
   block in code (element names in order of first appearance, coordinates
   in bohr). The file must not contain a `Geometry` block.
4. Writes the complete input to `dftbp_input.hsd`, then calls
   `setupCalculator`.
5. For each geometry: `setGeometry` → `getEnergy` (Mermin free energy, Ha)
   → `getGradients` (Ha/bohr).

## Output

| File | Content |
|---|---|
| screen | energies and gradients |
| `dftbp_result.dat` | energy and gradient at the `INxyz` geometry |
| `dftbp_input.hsd` | the full input DFTB+ received; renamed to `dftb_in.hsd` it can be run with the `dftb+` program for comparison |
| `dftbp.log` | DFTB+ screen output |

`IN_dftbplus.hsd` switches off the files DFTB+ normally writes
(`detailed.out`, `charges.bin`, `band.out`, `dftb_pin.hsd`).

**Required in `IN_dftbplus.hsd`:** `Analysis { PrintForces = Yes }`. In
DFTB+ 25.1 this is what makes the calculator compute forces when there is
no `Driver`; without it `getGradients` stops with "Forces not available".
(`Options { CalculateForces = Yes }` from DFTB+'s own API tests is the old
input format, ParserVersion 3, and is rejected with ParserVersion 14.)

## Results of step 3 (2026-11-17)

To repeat all four checks with one command (builds both programs, runs in a
new `tests/runs/verify-dftbplus_<timestamp>/`, prints PASS/FAIL):

```bash
tests/verify_dftbplus.sh
```

Pass criteria in the script: (a) |dE| < 1e-9 Ha, max |d grad| < 1e-10
Ha/bohr; (b) FD error < 1e-7 Ha/bohr at h = 1e-4 and at least 10× smaller
than at h = 1e-3; (c) |dE| < 1e-5 Ha, |d(energy change)| < 1e-6 Ha, max
|d grad| < 1e-4 Ha/bohr.

`h-c4h9-14` geometry (14 atoms, `INxyz` as is), `.spl` files from
`examples/sk/`, settings as in `IN_dftbplus.hsd`. Runs in
`tests/runs/dftbp-step3_20261117-181054/`.

| Check | How | Result |
|---|---|---|
| (d) `.spl` files accepted | driver run, `dftbp.log` | read without problems (incl. repeat counts like `19*0.0`); SCC converged in 8 iterations |
| (a) driver = `dftb+` program | `dftbp_input.hsd` run with `bin/dftb+` (only `WriteDetailedOut = Yes` changed), in `check_a/` | energy equal within 3e-11 Ha, gradient within 5e-13 Ha/bohr (printing precision); net force zero |
| (b) gradient = dE/dx | `dftbp_driver fd h`, in `check_b_h*/` | max \|analytic − FD\|: 2.3e-7 (h = 1e-3), 3.1e-9 (1e-4), 1.5e-9 (1e-5, noise floor) Ha/bohr; max gradient 1.5e-2 |
| (c) vs old engine | `tests/old_engine/old_engine` (calls `do_dftb` exactly as `derivs`), in `check_c/` | energy −10.708373138 (DFTB+) vs −10.708373482 (old): 3.4e-7 Ha; energy change for atom 1 +0.1 bohr: 1.3030151e-3 vs 1.3030623e-3 Ha; gradient max diff 1.2e-5, rms 2.8e-6 Ha/bohr |

Mermin free energy at the `INxyz` geometry: −10.708373137770 Ha. The
driver made about 170 DFTB+ calls in under 0.1 s in the finite-difference
runs.

## Status

- Step 1 (DFTB+ build): done.
- Step 2 (driver builds and links): done.
- Step 3 (checks a–d): done, all passed.
- Step 4 onward: see `agent_state.md`.
