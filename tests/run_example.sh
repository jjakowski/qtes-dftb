#!/bin/bash
# Run an example from ../examples in a fresh directory under tests/runs/.
#
#   usage:  tests/run_example.sh <example> [nprocs]
#   e.g.    tests/run_example.sh h-c4h9-14-ref 4
#
# The run directory gets the example's INdyn/INxyz/IN_dftb plus the shared
# Slater-Koster files from examples/sk/. Reference data is never touched.

set -e
here=$(cd "$(dirname "$0")" && pwd)
root=$(dirname "$here")

name=$1
np=${2:-4}
if [ -z "$name" ] || [ ! -d "$root/examples/$name" ] || [ "$name" = sk ]; then
  echo "usage: $0 <example> [nprocs]"
  echo "examples:" $(cd "$root/examples" && ls -d */ | tr -d / | grep -v '^sk$')
  exit 1
fi

exe=$root/qtes.x
[ -x "$exe" ] || { echo "missing $exe -- build it first (see examples/README.md)"; exit 1; }

run=$here/runs/${name}_$(date +%Y%m%d-%H%M%S)
mkdir -p "$run"
cp "$root/examples/$name"/{INdyn,INxyz,IN_dftb} "$run"/
cp "$root/examples/sk"/*.spl "$run"/
cp "$exe" "$run"/

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=1
ulimit -s unlimited 2>/dev/null || true

cd "$run"
echo "running $name with $np MPI ranks in $run"
mpirun -np "$np" ./qtes.x > stdout.log 2>&1
echo "done: $run/stdout.log"
