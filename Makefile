SHELL := /bin/bash
.SHELLFLAGS := -c

PYTHON_MODULE := python/3.14
VENV := .venv
ACTIVATE := source /etc/profile.d/lmod.sh && module load $(PYTHON_MODULE) && source $(VENV)/bin/activate

.PHONY: all env test clean

all: env

env: $(VENV)/bin/activate
	$(ACTIVATE) && pip install -e . --no-build-isolation && pip install pytest

$(VENV)/bin/activate: make_environment.sh requirements.txt
	bash make_environment.sh
	touch $@

test:
	python3 -m pytest tests/ -v

clean:
	rm -rf $(VENV) src/*.egg-info build dist
