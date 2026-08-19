PYTHON ?= python3
R_SCRIPT ?= Rscript

.PHONY: assemble check release-check r-setup figures draft clean

assemble:
	$(PYTHON) scripts/assemble.py

check:
	$(PYTHON) scripts/check.py

release-check:
	$(PYTHON) scripts/check.py --release

r-setup:
	R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null $(R_SCRIPT) --vanilla ../../benchmarks/figures/R/setup_packages.R

figures:
	R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null $(R_SCRIPT) --vanilla scripts/render_publication_figures_v3.R .
	R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null $(R_SCRIPT) --vanilla scripts/render_fastvep_cpu_complexity.R .
	$(PYTHON) scripts/record_provenance.py

draft: figures check
	$(PYTHON) scripts/build.py

clean:
	$(PYTHON) -c 'from pathlib import Path; import shutil; p=Path("build"); shutil.rmtree(p) if p.is_dir() else None'
