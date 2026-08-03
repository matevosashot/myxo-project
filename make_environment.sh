#!/usr/bin/env bash
set -euo pipefail

module load python/3.14

pip install --upgrade pip setuptools wheel
pip install -e . --no-build-isolation
