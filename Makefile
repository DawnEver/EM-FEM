.DEFAULT_GOAL := all

.PHONY: lint
lint:
	.venv/bin/ruff check src
	.venv/bin/ruff format --check src

.PHONY: fmt
fmt:
	.venv/bin/ruff check src --fix
	.venv/bin/ruff format src

.PHONY: test
test:
	coverage run -m unittest discover
	coverage html
	open htmlcov/index.html

.PHONY: all
all: lint test

.PHONY: clean
clean:
	rm -rf `find . -name __pycache__`
	rm -rf dist
	rm -rf build
	rm -rf .pytest_cache
	rm -rf .ruff_cache
	rm -rf .mypy_cache
	rm -rf htmlcov
	rm -rf *.egg-info
	rm -f .coverage
	rm -f .coverage.*
	rm -rf result.json
