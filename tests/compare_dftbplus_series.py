#!/usr/bin/env python3
"""Milestone 1, step 4: compare the series runs of tools/dftbplus/dftbp_driver
and tests/old_engine on the same geometries (a qtes-dftb qt.dat).

Usage: compare_dftbplus_series.py RUNDIR
RUNDIR holds the subfolders A_file, A_byqt, A_shuffle, fresh (driver runs:
one instance in three call orders, and a fresh instance per geometry) and
old (old engine). Missing subfolders are skipped.
"""
import os
import re
import sys


def read_driver(d):
    """Energies, times and call numbers in file order, gradients, SCC iterations."""
    rows = [l.split() for l in open(f"{d}/series_result.dat") if not l.startswith('#')]
    e = [float(r[4]) for r in rows]
    t = [float(r[5]) for r in rows]
    call = [int(r[3]) for r in rows]
    g = [list(map(float, l.split())) for l in open(f"{d}/series_grad.dat")]
    # SCC iterations per call, from the DFTB+ log (blocks appear in call order)
    log = open(f"{d}/dftbp.log").read().split('\n')
    iters, i = [], 0
    num = re.compile(r'^\s+\d+\s+-?\d\.\d+E[+-]\d+\s+-?\d\.\d+E[+-]\d+\s+\d\.\d+E[+-]\d+')
    while i < len(log):
        if 'iSCC' in log[i]:
            n = 0
            i += 1
            while i < len(log) and (num.match(log[i]) or log[i].startswith('>>')):
                if num.match(log[i]):
                    n += 1
                i += 1
            iters.append(n)
        else:
            i += 1
    scc = [iters[c - 1] if c - 1 < len(iters) else None for c in call]
    return dict(e=e, t=t, g=g, scc=scc)


def read_old(d):
    rows = [l.split() for l in open(f"{d}/old_series_result.dat") if not l.startswith('#')]
    g = [list(map(float, l.split())) for l in open(f"{d}/old_series_grad.dat")]
    return dict(e=[float(r[1]) for r in rows], t=[float(r[2]) for r in rows], g=g)


def diff(a, b):
    n = min(len(a['e']), len(b['e']))
    de = max(abs(a['e'][i] - b['e'][i]) for i in range(n))
    dg = max(abs(x - y) for i in range(n) for x, y in zip(a['g'][i], b['g'][i]))
    return n, de, dg


def mean(x):
    x = [v for v in x if v is not None]
    return sum(x) / len(x) if x else float('nan')


run = sys.argv[1]
R = {}
for name in ('A_file', 'A_byqt', 'A_shuffle', 'fresh'):
    if os.path.exists(f"{run}/{name}/series_result.dat"):
        R[name] = read_driver(f"{run}/{name}")
if os.path.exists(f"{run}/old/old_series_result.dat"):
    R['old'] = read_old(f"{run}/old")

print(f"Run folder: {run}")
print("\nPer call, mean over the series (for 'fresh' the time includes setting up the instance):")
print(f"  {'run':10s} {'geometries':>10s} {'SCC iter.':>10s} {'time (s)':>10s}")
for name, r in R.items():
    print(f"  {name:10s} {len(r['e']):10d} {mean(r.get('scc', [None])):10.2f} {mean(r['t']):10.4f}")

print("\nDifferences (max over all geometries):")
pairs = [('A_file', 'fresh', 'carry-over vs fresh start (1 rank, all QTs)'),
         ('A_byqt', 'fresh', 'carry-over vs fresh start (each QT in turn)'),
         ('A_shuffle', 'fresh', 'carry-over vs fresh start (shuffled order)'),
         ('A_file', 'A_shuffle', 'order dependence: file vs shuffled'),
         ('A_file', 'A_byqt', 'order dependence: file vs byqt'),
         ('fresh', 'old', 'DFTB+ (fresh) vs old engine'),
         ('A_file', 'old', 'DFTB+ (A, file order) vs old engine')]
for a, b, label in pairs:
    if a in R and b in R:
        n, de, dg = diff(R[a], R[b])
        print(f"  {label:45s} n={n:3d}  max|dE| = {de:.2e} Ha  max|d grad| = {dg:.2e} Ha/bohr")
if 'fresh' in R:
    gmax = max(abs(x) for g in R['fresh']['g'] for x in g)
    print(f"\nLargest gradient component in the series: {gmax:.2e} Ha/bohr")
