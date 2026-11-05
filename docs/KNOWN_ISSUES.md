# Known issues

This file lists open bugs and other behaviour worth knowing about. When an
issue is fixed, move it to **Fixed** and note the commit.

## Open

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

### 6. `ipot = 1` fit fails when too few trajectories remain

- **Where:** `qtes.f`, `subroutine qpot`. The linear quantum-force fit
  solves a least-squares system with `NBAS = nqd + 1` basis functions
  (4 for one quantum atom).
- **Problem:** the fit needs at least `NBAS` trajectories with nonzero
  weight, and in practice more, spread out in space. If too many
  trajectories are excluded below `qmn`, or if all trajectories coincide
  (`isr = 0`), the system is singular. The code then prints
  `DPOSV FAILED` and continues with a garbage quantum force, and the
  energy blows up.
- **Seen in:** a stress test (2026-11-05) of `h-flake-111` with `qmn`
  raised to 3.1 bohr. The surviving trajectories dropped 12 -> 8 -> 7 ->
  5 -> 3. At 3 remaining, the fit failed (t ~ 46), and the kinetic energy
  jumped to ~1e20 Ha. The pre-fix and post-fix binaries for issue F3
  behave the same way here.
- **Impact:** none of the examples are affected (`h-flake-111` uses the
  default `qmn = -15`).
- **Likely fix:** stop with a clear message when `DPOSV` fails, or when
  fewer than `NBAS` weighted trajectories remain, instead of continuing.
- **Found:** 2026-11-05

## Notes on behaviour (not bugs, but easy to miss)

- **`IEEE_DENORMAL` in restart runs is harmless.** gfortran reports at
  exit that some intermediate values were below about 2e-308
  (denormal). This comes from coordinates that are essentially zero
  (y about 1e-20 to 1e-26) being multiplied, and it does not affect the
  results.

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

### F3. `trot`: undefined weight index for trajectories below `qmn`

- **Where:** `qtes.f`, `subroutine trot`, force loop `do 24 ip=ip0,ip1`
- **Problem:** if a trajectory's quantum atom went below `qmn`, the
  loop set `w(n)=0d0`. `n` is never assigned in `trot` (it is an
  implicitly typed, uninitialized local), so the write went to an
  undefined position, possibly out of bounds, and the trajectory's own
  weight was not zeroed until the next half-step.
- **Fix:** `w(ip)=0d0` (2026-11-05), in the commit "Fix undefined weight
  index w(n) in trot (issue #1)". Note that `w` is not communicated
  between ranks after the force loop. The zero weight is therefore seen
  only on the owning rank until the next half-step zeroes it on all
  ranks.
- **Verification:**
  - All four examples are byte-identical before and after the fix,
    since none of them reaches this code path.
  - A stress test with `qmn = 3.1` bohr in `h-flake-111` exercises the
    path. Its results match the pre-fix binary up to the first
    in-step crossing (t = 18), then differ, as expected.
- **Remaining:** the momenta of excluded trajectories are still updated
  in the second half-step. This is harmless, because their weight is
  0.

### F4. Floating-point divide-by-zero (and invalid) at runtime

- **Symptom:** every run printed `IEEE_DIVIDE_BY_ZERO` at exit, and
  `h-flake-111` also printed `IEEE_INVALID_FLAG`.
- **Where:** `qtes.f`, `subroutine init_traj`, loop `do 91 i=1,nd`. The
  line `dx(i)=1.d0/dsqrt(4.d0*a0(i,i))` ran for all dimensions.
- **Problem:**
  - The wavepacket widths `a0` are set only for the quantum coordinates.
    For classical coordinates `a0(i,i)=0`, so `dx(i)` became
    `1/sqrt(0) = Inf`.
  - `a0` is dimensioned `a0(nqx,nqx)` with `nqx=200`. For systems with
    `nd > 200` (`h-flake-111`, `nd=333`), the loop read `a0(i,i)` out of
    bounds, and sqrt of a negative garbage value gave NaN.
  - `dx` is used only for `i <= nqd`, so the bad values never reached
    the results.
- **Found with:** a debug build (`-O0 -g -fbacktrace
  -ffpe-trap=zero,invalid`). After the fix, all four examples run with
  both traps active.
- **Fix:** compute `dx(i)` only inside the existing `if(i.le.nqd)`
  block (2026-11-05), in the commit "Fix divide-by-zero / out-of-bounds
  read in init_traj (issue #3)".
- **Verification:** all four examples are byte-identical before and
  after the fix, and the divide-by-zero and invalid messages are gone.
  The restart example still reports `IEEE_DENORMAL`, which is harmless
  (see the notes above).

### F5. `trot`: undefined forces for trajectories outside the z-window

- **Where:** `qtes.f`, `subroutine trot`, force loop `do 24 ip=ip0,ip1`
- **Problem:**
  - For a trajectory whose quantum atom was above `qmx` (or below
    `qmn`), the loop skipped DFTB with `goto 24`, which left
    `vlocal(ip)`/`dvlocal(:,ip)` unset.
  - These are uninitialized local arrays; with `-fopenmp` they live on
    the stack. They were still gathered and broadcast into `v`/`dv`, so
    the trajectory moved with undefined forces, which could be stale
    values or garbage.
  - With `ipf = 0`, the garbage also entered the Ehrenfest average,
    because the trajectory keeps its weight above `qmx`. It therefore
    corrupted the classical forces of every trajectory.
- **Fix:** outside the window, the code now sets the force to zero
  (`dvlocal(:,ip)=0`) and freezes the potential at its last value
  (`vlocal(ip)=v(ip)`; `v` holds the previous step on all ranks)
  (2026-11-05). This is free flight, as chosen by the code author.
  The below-`qmn` branch gets the same well-defined values, because its
  zero weight is visible only on the owning rank for one substep (F3).
- **Behaviour to know:**
  - The total energy jumps by a small step when a trajectory re-enters
    the window and DFTB resumes. The size of the step is how much the
    true potential changed while it was frozen.
  - Choose `qmx` where the potential is flat. In a deliberately harsh
    test, `h-c4h9-14-ref` with `qmx = 8.5` bohr (about 6 bohr above the
    C-H group), the step was -2.2e-5 Ha. The H flew at constant momentum
    until z < 8.5, and then followed DFTB forces.
- **Verification:** all four examples are byte-identical before and
  after the fix, since none of them leaves the window.
