# Dynamics modes and switches

This file describes what the code actually does for each run-mode switch
in `INdyn` and `INxyz`, and where in the source it does it. Line numbers
refer to the input files. For bugs in these code paths, see
`KNOWN_ISSUES.md`.

## Overview

- Each **quantum trajectory** (QT) `i` is a full copy of the system:
  quantum atoms, then classical atoms. It has positions `q(:,i)`, momenta
  `p(:,i)`, an action `s(i)`, and a weight `w(i)`.
- The quantum atoms are listed first in `INdyn`. All current examples
  have one quantum atom, the incoming H (atom 1).
- `trot` (`qtes.f`) propagates every trajectory with velocity Verlet.
  Each large step (`dt`, `INdyn` line 5) is split into `kout` substeps
  (line 4), and there are `kmax` large steps (line 3).
- At every substep, each trajectory gets one DFTB calculation for its own
  geometry, which gives its potential `v(i)` and forces `dv(:,i)`.
  Trajectories are distributed over the MPI ranks.
- Units are atomic units: bohr, hartree, and electron masses (1836 for H).
  `INxyz` gives the geometry in Angstrom, which the code converts to bohr
  on input.

## Switches

| Switch | Input | Values | Code |
|---|---|---|---|
| `ipot` | `INdyn` line 6, 1st value | 0 = no quantum potential; 1 = linear quantum-force fit in all quantum dims; 2 = linear fit along z only | `qpot`, `trot` |
| `isr` | `INdyn` line 6, 2nd value | 0 = no sampling (all QTs at the wavepacket centre); 1 = Gaussian sampling | `init_traj` |
| `pow` | `INdyn` line 7 | sampling cutoff: keep points with Gaussian weight > 10^-pow | `init_traj` |
| `ipf` | `INdyn` line 8 | 0 = Ehrenfest; any other value = independent classical copies | `init_traj`, `trot`, `set_traj` |
| `np`, restart | `INdyn` line 9, values 1-2 | number of QTs; 0 = start fresh, 1 = restart from `traj.sav` | `read_input_mpi`, `set_traj` |
| `nqtprt`, `nclprt` | `INdyn` line 9, optional values 3-4 | QTs printed in `qtx..qtz` / classical atoms in `clx..clz` (default 20 each; 0 or less = all) | `read_input_mpi`, `aver` |
| `qmn`, `qmx` | `INxyz` line 6 | z-window (bohr) for the quantum atom | `trot` |

All current examples use `ipot=0`, `isr=0` and `ipf=0`, except
`h-flake-111`, which uses `ipot=1` and `isr=1` (a chemisorbed H with the
quantum potential and sampling).

## `ipf`: Ehrenfest vs independent classical copies

`ipf` controls one thing only: whether the forces on the **classical**
atoms are averaged over the trajectories. The averaging block appears in
three places: `init_traj` (initial forces), `trot` (every substep), and
`set_traj` (restart). The forces on the quantum atoms are never averaged.

### `ipf = 0`: Ehrenfest mode

- The stdout log prints `EHRENFEST: AVERAGE OVER QDOF 0`.
- Every classical force is replaced by the weighted average over all
  trajectories, `<F> = sum_i w_i F_i / sum_i w_i`, and all trajectories
  get the same value.
- All trajectories start with the same classical positions and zero
  classical momenta. Together with the shared force, this keeps the
  **classical atoms identical in every trajectory**. They move as one
  mean-field system that responds to the averaged quantum atom.
- Each trajectory's quantum atom moves under its own DFTB force, plus the
  quantum force when `ipot` is not 0.

### `ipf != 0`: independent classical copies

- The stdout log prints `DO NOT AVERAGE OVER QDOF <ipf>`. All nonzero
  values behave the same.
- No averaging is done. Each trajectory carries its own copy of the
  classical atoms, which move under the DFTB force for that trajectory's
  own geometry. **Classical atoms therefore differ between
  trajectories.** Where the H is in trajectory *i* determines how the
  target moves in trajectory *i*.
- The trajectories are coupled only through the quantum potential
  (`ipot != 0`), and only in the quantum coordinates. With `ipot = 0`,
  the trajectories are fully independent.

### Effect on output

- In both modes, `qtx/qty/qtz` give the quantum atom's x/y/z coordinate
  for up to `nqtprt` QTs.
- What `clx/cly/clz` contain depends on the mode (`aver` in `qtes.f`):
  - **`ipf = 0`:** the classical atoms are the same in every QT, so the
    files give up to `nclprt` classical atoms, one atom per column. The
    values are taken from the first QT with nonzero weight.
  - **`ipf != 0`:** the classical atoms differ between QTs, so the files
    give the first classical atom, one QT per column, for up to `nqtprt`
    QTs.

### Caveats for both modes

- With `isr = 0`, all trajectories start at the same point. They then
  stay identical even when `ipf != 0`, so the choice of `ipf` only
  matters when `isr = 1`.
- A trajectory excluded below `qmn` stops moving, classical atoms
  included. See the `qmn` / `qmx` section below.
- No example uses `ipf != 0`, so this mode has not been tested with the
  current build.

## `ipot`: quantum potential

- **`ipot = 0`:** `qpot` returns zero, so there is no quantum
  potential or force. The quantum atom moves classically under its own
  DFTB force.
- **`ipot = 1`:** the quantum force comes from a linear fit of the
  momentum, across the ensemble, in all quantum coordinates. The fit
  coefficients are written to `pfit`.
- **`ipot = 2`:** the same linear fit, but only along z. In `trot`, the
  x and y components of the quantum atom get only the DFTB force, and z
  gets the DFTB force plus the quantum force.
- With `isr = 0`, all trajectories are identical, so the fitted quantum
  force carries no information. The quantum potential only makes sense
  with sampling on.
- Likely bug: with `ipot = 1`, the force on the quantum atom is
  applied twice in each half-step (`KNOWN_ISSUES.md` #5).

## `isr`, `pow`: initial sampling

- **`isr = 1`:** each quantum coordinate is sampled as
  `q = q0 + dx * gaussian_deviate`, with `dx = 1/sqrt(4 a0)`, where `a0`
  is the width parameter in `INdyn`.
  - Points whose Gaussian weight is at or below `10^-pow` are rejected.
  - Each accepted point gets weight `w = 1/np`.
  - All momenta are set to `p0`, with no momentum sampling.
- **`isr = 0`:** every trajectory is placed at the centre, `q0` and `p0`.
- Classical atoms are never sampled. They always start at the `INxyz`
  geometry plus the `INdyn` displacement, with **zero momentum**. The
  classical px/py/pz values in `INdyn` are ignored.

## `qmn`, `qmx`: z-window for the quantum atom

`qmn` and `qmx` do not act on individual atoms. In every substep, `trot`
checks the z-coordinate of the **quantum atom (atom 1)** in each
trajectory, `q(3,i)`, and applies the result to that **whole
trajectory**:

| Quantum-atom z in trajectory *i* | What `trot` does | Intended meaning |
|---|---|---|
| `qmn <= z <= qmx` | normal propagation | the atom is in the interaction region |
| `z < qmn` | sets `w(i) = 0`, stops updating all positions in the trajectory, and skips DFTB | the atom has passed through or beyond the target, so the trajectory is **removed from the ensemble**. It no longer contributes to averages, energies, or the Ehrenfest force, and the norm in `aver` (column 2) drops below 1 |
| `z > qmx` | skips DFTB; the force is set to zero and the potential is frozen at its last value. The trajectory keeps moving with its weight (free flight) | the atom is far above the target, so the expensive DFTB call is avoided |

Notes:

- An equivalent `z > qmx` check in the position update is commented out
  (`qtes.f`, first half-step of `trot`).
- Free flight above `qmx` means the total energy takes a small step
  when a trajectory re-enters the window and DFTB resumes. Choose `qmx`
  where the potential is flat (see `KNOWN_ISSUES.md` F5).
- Known issues in this path:
  - #6: with `ipot = 1`, the quantum-force fit fails if too few
    trajectories remain inside the window.
  - F3 (fixed): the force loop zeroed an undefined weight `w(n)`.
    It now zeroes `w(ip)`.
  - F5 (fixed): forces outside the window were undefined. They are now
    zero, with the potential frozen.
  - Momenta of excluded trajectories still change in the second
    half-step. This is harmless, because their weight is 0.
- Most examples use `-15 15` bohr; `h-flake-53` uses `-15 50`. In the
  current examples the H stays well inside this window. For example, in
  `h-c4h9-14-ref` it moves from z = 8.8 to 8.0 bohr.

## Restart

A restart is selected with the second value on `INdyn` line 9 set to 1.

- `read_input_mpi` reads the time `t` and the number of trajectories from
  `traj.sav`, and shortens the run to the remaining `kmax - t/dt` large
  steps.
- `set_traj` reloads all positions, momenta, actions and weights, and
  then recomputes the forces, including the Ehrenfest average when
  `ipf = 0`.
- The output files are opened with `status='old'` and appended to, so
  they must already exist.
- `h-c4h9-14-restart` tests this path. It reproduces a single
  uninterrupted run to round-off.
