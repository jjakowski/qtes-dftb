#!/bin/bash
# Milestone 1, step 3: verify that tools/dftbplus/dftbp_driver gets correct
# energies and forces from DFTB+ (checks a-d, see tools/dftbplus/README.md).
#
# Usage:  tests/verify_dftbplus.sh
# Builds the driver and the old-engine test program, runs everything in a
# new folder tests/runs/verify-dftbplus_<timestamp>/ and prints PASS/FAIL
# per check. Exit code 0 only if all checks pass.
#
# DFTB+ install: DFTBPLUS_PREFIX (default as in tools/dftbplus/Makefile).

REPO=$(cd "$(dirname "$0")/.." && pwd)
DFTBPLUS_PREFIX=${DFTBPLUS_PREFIX:-/home/j2c/Development/DFTB/dftbplus-25.1/install-qtes}
EX=$REPO/examples/h-c4h9-14
SK=$REPO/examples/sk
SPL="cc.spl ch-g2ms.spl hc-g2ms.spl hh.spl"     # the four files IN_dftbplus.hsd names

source ~/miniconda3/etc/profile.d/conda.sh && conda activate rmg-dev
set -u            # after conda activate: its scripts use unset variables
ulimit -s unlimited
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

echo "Building driver and old-engine test program ..."
make -s -C "$REPO/tools/dftbplus" DFTBPLUS_PREFIX="$DFTBPLUS_PREFIX" > /dev/null || exit 1
make -s -C "$REPO/tests/old_engine" > /dev/null || exit 1
DRV=$REPO/tools/dftbplus/dftbp_driver
OLD=$REPO/tests/old_engine/old_engine

V=$REPO/tests/runs/verify-dftbplus_$(date +%Y%m%d-%H%M%S)
mkdir -p "$V" && cd "$V" || exit 1
echo "Run folder: $V"

# ---- driver run (also check d): INxyz, INdyn, IN_dftbplus.hsd, 4 .spl files
mkdir drv && cp "$EX/INxyz" "$EX/INdyn" "$REPO/tools/dftbplus/IN_dftbplus.hsd" drv/
for f in $SPL; do cp "$SK/$f" drv/; done
(cd drv && "$DRV" > driver.out 2>&1); RC_DRV=$?

# ---- (a) the same input through the dftb+ program
mkdir check_a
for f in $SPL; do cp "$SK/$f" check_a/; done
sed 's/WriteDetailedOut = No/WriteDetailedOut = Yes/' drv/dftbp_input.hsd > check_a/dftb_in.hsd
(cd check_a && "$DFTBPLUS_PREFIX/bin/dftb+" > dftb.out 2>&1); RC_A=$?

# ---- (b) finite differences, two step sizes
for h in 1e-3 1e-4; do
  mkdir check_b_h$h && cp drv/INxyz drv/INdyn drv/IN_dftbplus.hsd check_b_h$h/
  for f in $SPL; do cp "$SK/$f" check_b_h$h/; done
  (cd check_b_h$h && "$DRV" fd $h > driver_fd.out 2>&1)
done

# ---- (c) old engine on the same geometry (also needs IN_dftb)
mkdir check_c && cp "$EX/INxyz" "$EX/INdyn" "$EX/IN_dftb" check_c/
for f in $SPL; do cp "$SK/$f" check_c/; done
(cd check_c && mpirun -np 1 "$OLD" > old.out 2>&1); RC_C=$?

# ---- evaluate
RC_DRV=$RC_DRV RC_A=$RC_A RC_C=$RC_C python3 - <<'EOF'
import os, re, sys

def result(path):
    lines = [l for l in open(path) if l.strip() and not l.startswith('#')]
    return float(lines[0]), [list(map(float, l.split()[1:4])) for l in lines[1:]]

def energies(path):
    return [float(x) for x in re.findall(r'Mermin free energy:\s+(\S+)', open(path).read())]

ok_all = True
def check(name, ok, detail):
    global ok_all
    ok_all &= ok
    print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")

# (d) driver ran, DFTB+ read the .spl files
log = open('drv/dftbp.log').read() if os.path.exists('drv/dftbp.log') else ''
check("(d) .spl files read, driver ran", os.environ['RC_DRV'] == '0' and 'Reading SK-files' in log,
      f"driver exit {os.environ['RC_DRV']}")
if os.environ['RC_DRV'] != '0':
    sys.exit(1)
E, g = result('drv/dftbp_result.dat')
nat = len(g)

# (a) driver = dftb+ program
d = open('check_a/detailed.out').read() if os.environ['RC_A'] == '0' else ''
if d:
    Ea = float(re.search(r'Total Mermin free energy:\s+(\S+) H', d).group(1))
    i = d.index('Total Forces')
    F = [list(map(float, l.split()[1:4])) for l in d[i:].split('\n')[1:nat + 1]]
    dE = abs(E - Ea)
    dg = max(abs(g[a][k] + F[a][k]) for a in range(nat) for k in range(3))
    check("(a) driver = dftb+ program", dE < 1e-9 and dg < 1e-10,
          f"|dE| = {dE:.1e} Ha (< 1e-9), max |d grad| = {dg:.1e} Ha/bohr (< 1e-10)")
else:
    check("(a) driver = dftb+ program", False, f"dftb+ exit {os.environ['RC_A']}")

# (b) finite differences
err = {}
for h in ('1e-3', '1e-4'):
    t = open(f'check_b_h{h}/driver_fd.out').read()
    m = re.search(r'max \|analytic - FD\| =\s+(\S+)', t)
    err[h] = float(m.group(1)) if m else float('inf')
check("(b) gradient = dE/dx", err['1e-4'] < 1e-7 and err['1e-4'] < err['1e-3'] / 10,
      f"max |analytic - FD| = {err['1e-3']:.1e} (h=1e-3), {err['1e-4']:.1e} (h=1e-4) Ha/bohr; "
      f"needs < 1e-7 and > 10x smaller with the smaller step")

# (c) old engine
if os.environ['RC_C'] == '0':
    Eo, go = result('check_c/old_result.dat')
    en, eo = energies('drv/driver.out'), energies('check_c/old.out')
    dE = abs(E - Eo)
    dDE = abs((en[1] - en[0]) - (eo[1] - eo[0]))
    dg = max(abs(g[a][k] - go[a][k]) for a in range(nat) for k in range(3))
    check("(c) DFTB+ vs old engine", dE < 1e-5 and dDE < 1e-6 and dg < 1e-4,
          f"|dE| = {dE:.1e} Ha (< 1e-5), |d(energy change)| = {dDE:.1e} Ha (< 1e-6), "
          f"max |d grad| = {dg:.1e} Ha/bohr (< 1e-4)")
else:
    check("(c) DFTB+ vs old engine", False, f"old engine exit {os.environ['RC_C']}")

print("ALL CHECKS PASSED" if ok_all else "SOME CHECKS FAILED")
sys.exit(0 if ok_all else 1)
EOF
