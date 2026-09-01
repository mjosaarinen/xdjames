# D-James / James -- see README.md, and ref/README.md for the upstream material.

PY    ?= python3
SAGE  ?= sage
CARGO ?= cargo
PKG   := djames-py
RS    := djames-rs

.PHONY: help test test-full kat fieldpolys demo clean \
        rs-test rs-test-full rs-bench rs-fmt rs-clippy check

help:
	@echo "targets:"
	@echo "  test        test suite at toy parameters, re-deriving kat/toy.json"
	@echo "  test-full   as above, also re-deriving kat/q2.json and kat/keygen.json"
	@echo "  kat         regenerate every known-answer test vector"
	@echo "  fieldpolys  regenerate djames/fieldpoly.json from scratch"
	@echo "  demo        run the authors' SageMath proof-of-concept (needs sage)"
	@echo "  clean       remove Python, Sage and Rust build artifacts"
	@echo ""
	@echo "  rs-test      Rust unit tests + toy KAT vectors"
	@echo "  rs-test-full Rust tests over all 20 real parameter sets"
	@echo "  rs-bench     Rust keygen/sign/verify timings"
	@echo "  rs-fmt       rustfmt"
	@echo "  rs-clippy    clippy, warnings denied"
	@echo "  check        both implementations, on the same vectors"

test:
	cd $(PKG) && $(PY) -m unittest discover -s tests -v

test-full:
	cd $(PKG) && DJAMES_FULL_KAT=1 $(PY) -m unittest discover -s tests -v

kat:
	cd $(PKG) && $(PY) tools/gen_kat.py toy
	cd $(PKG) && $(PY) tools/gen_kat.py q2
	cd $(PKG) && $(PY) tools/gen_kat.py keygen

fieldpolys:
	cd $(PKG) && $(PY) tools/gen_fieldpolys.py

demo:
	cd ref && $(SAGE) djames_demo.sage

# --- Rust ------------------------------------------------------------------
# The Rust crate reads djames-py/kat/*.rsp, so `make kat` must have run first
# if the vectors ever change.

rs-test:
	cd $(RS) && $(CARGO) test --release

rs-test-full:
	cd $(RS) && $(CARGO) test --release -- --include-ignored

rs-bench:
	cd $(RS) && $(CARGO) run --release --example bench

rs-fmt:
	cd $(RS) && $(CARGO) fmt

rs-clippy:
	cd $(RS) && $(CARGO) clippy --release --all-targets -- -D warnings

# Both implementations against the same known-answer vectors.
check: test rs-test

# Removes only regenerable scratch.  The checked-in artifacts -- kat/*.json
# and djames/fieldpoly.json -- are left alone; use the targets above to
# rebuild those deliberately.
clean:
	@echo "cleaning Python artifacts"
	@find . -name '__pycache__' -prune -exec rm -rf {} +
	@find . -name '*.py[cod]' -delete
	@for n in .pytest_cache .mypy_cache .ruff_cache .tox .coverage htmlcov \
	          build dist '*.egg-info'; do \
		find . -name "$$n" -prune -exec rm -rf {} + ; \
	done
	@echo "cleaning Sage artifacts"
	@find . -name '*.sage.py' -delete
	@find . -name '*.sobj' -delete
	@find . -name '*.spyx.c' -delete
	@find . -name '*.spyx.so' -delete
	@find . -name '.sage' -type d -prune -exec rm -rf {} +
	@echo "cleaning Rust artifacts"
	@find . -name 'target' -type d -prune -exec rm -rf {} +
	@find . -name 'Cargo.lock' -delete
	@echo "cleaning notebook checkpoints"
	@find . -name '.ipynb_checkpoints' -type d -prune -exec rm -rf {} +
	@echo "done"
