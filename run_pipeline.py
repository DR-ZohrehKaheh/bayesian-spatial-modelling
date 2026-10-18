#!/usr/bin/env python3
# Project supervisor: Dr. Zohreh Kaheh.
# Undergraduate project collaborator: Jingxuan Hu.
"""Compile, run and verify the native simulation workflow."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('outputs'))
    parser.add_argument('--compiler', default='g++')
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    out = args.output.resolve()
    build = out / 'build'
    build.mkdir(parents=True, exist_ok=True)
    if not shutil.which(args.compiler):
        raise RuntimeError('A C++17 compiler is required')
    for name in ['test_bd_secr_core', 'run_multidataset_recovery']:
        exe = build / name
        subprocess.run([args.compiler, '-O3', '-std=c++17', '-I'+str(source),
                        str(source / (name+'.cpp')), '-o', str(exe)], check=True)
        command = [str(exe)]
        if name == 'run_multidataset_recovery':
            command.append(str(out / 'results/full_recovery'))
        subprocess.run(command, check=True)
    subprocess.run([sys.executable, str(source/'summarize_recovery.py'), str(out)], check=True)
    results = out / 'results/full_recovery'
    traces = pd.read_csv(results / 'recovery_traces.csv')
    acceptance = pd.read_csv(results / 'recovery_acceptance.csv')
    groups = traces.groupby(['dataset', 'chain']).size()
    if len(traces) != 72000 or len(groups) != 12 or not groups.eq(6000).all():
        raise RuntimeError('Incomplete simulation output')
    if not np.isfinite(traces.select_dtypes('number').to_numpy()).all():
        raise RuntimeError('Non-finite simulation output')
    if not (traces.N.ge(0).all() and traces.sigma.gt(0).all()):
        raise RuntimeError('Invalid parameter support')
    if not acceptance.cache_max_abs_error.lt(1e-8).all():
        raise RuntimeError('Intensity cache check failed')
    print('Native pipeline execution and output-integrity checks passed.')


if __name__ == '__main__':
    main()

