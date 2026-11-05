# Known issues

This file lists open bugs and other behaviour worth knowing about. When an
issue is fixed, move it to **Fixed** and note the commit.

## Open

### 1. `trot`: wrong weight zeroed for trajectories below `qmn`

- **Where:** `qtes.f`, `subroutine trot`, force loop `do 24 ip=ip0,ip1`
- **Problem:** if a trajectory's quantum atom goes below `qmn` (the
  z-range set in `INxyz`), the loop runs `w(n)=0d0` and skips that
  trajectory. The loop index is `ip`, and `n` has an unrelated value at
  that point, so the wrong weight gets zeroed, or none at all. The
  half-step loop earlier in `trot` uses the correct index (`w(i)=0d0`).
- **Impact:** this only happens when a trajectory goes below `qmn`.
  None of the current examples reach that point, so their results are
  not affected. Also, the momenta of excluded trajectories are still
  updated in the second half-step. That is harmless, because their
  weight is 0.
- **Likely fix:** `w(ip)=0d0`. This is the same kind of bug as fixed
  issue F1.
- **Found:** 2026-11-05

### 3. Floating-point divide-by-zero at runtime

- **Symptom:** every run prints `IEEE_DIVIDE_BY_ZERO` at exit, in all
  examples.
- **Impact:** none observed. The results are finite and match the
  references.
- **To investigate:** build with `-ffpe-trap=zero -g` to find the line
  where it happens.

### 4. `trot`: undefined forces for trajectories above `qmx`

- **Where:** `qtes.f`, `subroutine trot`, force loop:
  `if(q(3,ip).gt.qmx) goto 24`
- **Problem:** for a trajectory whose quantum atom is above `qmx`, the
  loop skips the DFTB call. `vlocal(ip)` and `dvlocal(:,ip)` are then
  left unset. These are local arrays that are never initialized, and
  with `-fopenmp` they live on the stack. The code still gathers them
  with `mpi_gather` and uses them as that trajectory's potential and
  forces, so the trajectory gets garbage values.
- **Impact:** this only happens when a trajectory's quantum atom goes
  above `qmx`. None of the current examples reach that point.
- **Likely fix:** this needs the author's intent. One option is to set
  the potential and force to zero (free flight far from the target).
  Another is to keep the previous step's values, for example by copying
  `v`/`dv` into `vlocal`/`dvlocal` before the loop.
- **Found:** 2026-11-05

### 5. `trot`: with `ipot = 1`, the quantum-atom force is applied twice (likely bug)

- **Where:** `qtes.f`, `subroutine trot`, both momentum half-steps
  (the lines `if(ipot.eq.1) p(k,i)=p(k,i)-(dv(k,i)+du(k,i))*dt2`)
- **Problem:** for each quantum coordinate, the `ipot.eq.1` line
  subtracts `(dv+du)*dt2`. The loop over components that follows then
  subtracts the force again: `dv*dt2` for x and y, and `(dv+du)*dt2`
  for z. The total kick is `2dv+du` in x and y and `2(dv+du)` in z,
  instead of `dv+du`. The per-component loop was probably meant for
  `ipot != 1` only, as an `else` branch.
- **Impact:** this affects only `ipot = 1` runs. Of the examples, only
  `h-flake-111` uses `ipot = 1`. Its results, including its reference
  data if one is added, may therefore contain this error.
- **Status:** likely bug. The code author agrees (2026-11-05). The fix is
  deferred, and it blocks the `h-flake-111` reference.
- **Found:** 2026-11-05

## Notes on behaviour (not bugs, but easy to miss)

See `docs/DYNAMICS.md` for a full description of the run modes (`ipf`,
`ipot`, `isr`, `qmn`/`qmx`, restart).

- **Classical initial momenta are ignored.** `init_traj` sets the
  initial momenta of classical atoms to 0, so the px/py/pz values that
  `INdyn` gives for classical atoms are never used. All current examples
  use 0 for them anyway.
- **With `isr = 0`, all trajectories are identical.** Sampling is off,
  so every trajectory starts at the wavepacket centre. All current
  examples use `isr = 0`.
- **With Ehrenfest averaging (`ipf = 0`), classical atoms are identical
  across trajectories.** Classical forces are averaged over the
  trajectories in `trot`. The exception is a trajectory excluded below
  `qmn`, which stops moving.
- **Restart needs the old output files.** A restart (flag 1 on line 9
  of `INdyn`) opens the output files with `status='old'` and appends to
  them, so they must already exist (see `examples/h-c4h9-14-restart`).

## Fixed

### F1. `aver`: `qtx/qty/qtz` and `clx/cly/clz` contained only the time column

- **Where:** `qtes.f`, `subroutine aver`
- **Problem:** the trajectory count is passed in as `n`, but the code
  used `np`, which is undeclared and uninitialized, in
  `ip = min(np,20)`. As a result, no trajectory columns were printed.
- **Fix:** `ip = min(n,20)`, in commit `ccf8435` (2026-11-05). The
  reference outputs were regenerated.

### F2. `h-flake-111` example could not run

- **Where:** `examples/h-flake-111/INdyn`
- **Problem 1:** line 6 had only `ipot` (`1`), but the code reads two
  values, `ipot` and `isr`, so the run stopped with `Bad integer for
  item 2 in list input`. The file predates the `isr` flag.
- **Problem 2:** the H line shifted the H by z0 = -2.0 bohr. `INxyz`
  already places this H chemisorbed, 1.118 A above C2, so the shift
  left it 0.06 A from C2. DFTB then failed on the first step
  (`ewevge: ier = 746`, overlap matrix not positive definite). The
  `INdyn` file was probably written for a different `INxyz` geometry.
- **Problem 3:** with `isr = 0`, all trajectories are identical, so
  the `ipot = 1` quantum-force fit is singular (`DPOSV FAILED`).
- **Fix:** line 6 changed to `1 1` (`ipot=1`, sampling on), and the H
  shift changed to z0 = 0. The example now simulates a chemisorbed H
  with quantum effects. The author confirmed this was the intent
  (2026-11-05). Fixed in the commit "Fix h-flake-111 inputs: chemisorbed
  H with sampling on".
