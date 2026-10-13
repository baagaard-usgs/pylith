#!/usr/bin/env python3

import subprocess

commands = [
    "./generate_gmsh.py --write",
    "pylith step01_inflation.cfg",
    "pylith step01_inflation.cfg step01b_inflation.cfg",
    "pylith step02_inflation_statevars.cfg",
]

for cmd in commands:
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR ({result.returncode}): {cmd}")
        print(result.stderr.strip())
    else:
        print(f"OK: {cmd}")
