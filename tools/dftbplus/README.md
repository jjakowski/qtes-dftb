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
(e.g. from `examples/sk/`):

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

## Status

- Step 1 (DFTB+ build): done.
- Step 2 (driver builds and links): done.
- Step 3 onward (checks): see `agent_state.md`.
