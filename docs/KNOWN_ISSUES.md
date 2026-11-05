# Known issues

This file lists open bugs and other behaviour worth knowing about. When an
issue is fixed, move it to **Fixed** and note the commit.

## Open

No open issues (2026-11-05).

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
  index w(n) in trot (issue #1)". That zero weight was still seen only
  on the owning rank. This was fixed later, in F7.
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
  The below-`qmn` branch gets the same well-defined values. Since F7 its
  weight is already zero on every rank, so the values are only a
  safeguard.
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

### F6. Quantum-force fit failure ended the run with a rank-local `stop`

- **Where:** `qtes.f`, `subroutine qpot` (`ipot = 1` and `ipot = 2`)
- **What the fit is:**
  - The linear quantum-force (LQF) method fits the non-classical
    momentum with the basis `f = (x, y, z, 1)`, i.e. `NBAS = nqd + 1`
    functions (for `ipot = 2`, only z and 1).
  - It builds the weighted moment matrix `M = sum_k w_k f(q_k) f(q_k)^T`
    over all QTs and solves it with LAPACK `dposv`.
  - `M` can be inverted only if the QTs with `w > 0` span all fitted
    directions: at least `NBAS` of them, not lying on a common plane,
    line, or point.
- **Problem:**
  - When `M` was singular, `qpot` printed `DPOSV FAILED` and called
    `stop` on that rank only. The other ranks hung until `mpirun`
    killed the job with a generic error.
  - `M` is singular when fewer than `NBAS` QTs remain inside the
    `qmn`/`qmx` window, or when all QTs coincide (`ipot != 0` with
    `isr = 0`).
- **Fix** (2026-11-05), in the commit "Abort cleanly when the
  quantum-force fit cannot be done (issue #6)":
  - `qpot` now counts the QTs with `w > 0` before the solve and stops if
    there are fewer than `NBAS`.
  - If `dposv` fails, `qpot` stops as well.
  - In both cases the new subroutine `qpfail` prints a clear
    `QPOT ERROR` (time, count, `info`, likely cause, hints) and calls
    `mpi_abort`, so all ranks stop together.
- **Remaining limitation (not a bug):** before `M` becomes exactly
  singular, a fit from few QTs, or from QTs that nearly lie on a common
  plane or line, is ill-conditioned. The quantum force can then become
  very large.
  - In the `qmn = 3.1` stress test, the quantum kinetic energy reached
    about 11.6 Ha with 5 QTs left, before the count check stopped the
    run at 3 QTs (t = 46.5).
  - Use enough QTs, with `isr = 1`, and a `qmn`/`qmx` window that keeps
    most of them inside.
- **Verification:**
  - All four examples are byte-identical before and after the fix.
  - In the stress test, all 12 ranks report `only 3 QTs with w>0; fit
    needs 4` at t = 46.5, followed by one clean `MPI_ABORT`.
  - The run with `ipot = 1` and `isr = 0` reports `DPOSV FAILED, info=3
    ... moment matrix singular ... (isr=0?)` at t = 0.

### F7. Weights differed between MPI ranks after a QT crossed `qmn`

- **Where:** `qtes.f`, `subroutine trot`.
- **Background:** every rank keeps its own full copy of `q`, `p`, `s`
  and `w` for all QTs. The copies stay identical only because all ranks
  make the same updates. They are sent between ranks once, at the
  start, and never again.
- **Problem:**
  - When a QT's quantum atom crossed below `qmn` during a substep, only
    the rank that owns that QT set `w(ip)=0`, in the force loop.
  - The Ehrenfest average that follows uses each rank's own `w`. The
    owning rank therefore averaged the classical forces without that QT,
    while every other rank still included it.
  - As a result, the classical-atom momenta of all QTs were updated
    differently on different ranks. From then on the ranks' copies of
    the system drifted apart for the rest of the run. The next half-step
    made `w` agree again, but `q` and `p` were never re-synchronized.
- **Fix:** in the first half-step loop, which every rank runs for all
  QTs, `w(i)=0` is now also set right after the move when
  `q(3,i) < qmn` (2026-11-05), in the commit "Exclude QTs crossing qmn
  on every rank (issue #7)". All ranks make this decision from
  identical data, so no communication is needed. A crossing QT is
  dropped from `qpot` and from the Ehrenfest average in the same
  substep in which it crosses. The `w(ip)=0` in the force loop stays as
  a harmless safeguard.
- **Verification:**
  - Scratch builds wrote a per-rank checksum (`sum q`, `sum p`,
    `sum w`) after every substep. The test case was `h-flake-111` with
    `qmn = 3.1` bohr on 12 ranks.
  - Old code: the ranks were identical until the first in-step crossing
    at t = 18. There `sum w` took 2 different values, and `sum q` and
    `sum p` differed by up to 4e-4 relative, a difference that persisted
    until the end of the run.
  - Fixed code: all 12 ranks stayed bit-identical at every substep.
  - All four examples are byte-identical before and after the fix,
    since none of them crosses `qmn`.

### F8. `trot`: with `ipot = 1`, the quantum-atom force was applied twice

- **Where:** `qtes.f`, `subroutine trot`, both momentum half-steps.
- **Problem:**
  - For each quantum coordinate, the line
    `if(ipot.eq.1) p=p-(dv+du)*dt2` gave a full kick.
  - The per-component loop that follows ran for every `ipot`, and it
    kicked the same coordinate again: `dv` in x and y, `dv+du` in z.
  - With `ipot = 1` the total kick was therefore `2dv+du` in x and y,
    and `2(dv+du)` in z. The H felt twice the DFTB force, while the
    energy written to `ener` assumes a single force, so the equations of
    motion did not match the reported energy.
  - The loop alone is correct for `ipot = 0` and `ipot = 2` (QP along
    z). The `ipot = 1` line had been added without making the two
    alternatives.
- **Fix:** the `ipot = 1` line and the per-component loop are now
  `if/else` alternatives (2026-11-05), in the commit "Apply the
  quantum-atom force once with ipot=1 (issue #5)". Every quantum
  coordinate now gets exactly one kick per half-step. `ipot = 0` and
  `ipot = 2` are unchanged.
- **Evidence:** energy conservation, max |E - E0|, before -> after the
  fix:

  | Case | Double kick | Single kick |
  |---|---|---|
  | `h-flake-111` (chemisorbed H, 12 QTs, t = 0 to 75) | 1.3e-3 Ha | < 1e-4 Ha (no change in the 7 printed digits) |
  | `h-c4h9-14-ref` inputs with `ipot=1 isr=1` (incoming H, 48 QTs, t = 0 to 100) | 2.5e-3 Ha, steady drift | 2.3e-6 Ha |

- **Verification:**
  - The three `ipot = 0` examples are byte-identical before and after
    the fix, with their references unchanged.
  - `h-flake-111` changed as intended.
  - Its first reference was created from the fixed code.
