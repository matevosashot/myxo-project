# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment setup

Development happens on an HPC cluster. Activate the venv with:

```bash
module load python/3.14
source .venv/bin/activate

pip install -e . --no-build-isolation

```
