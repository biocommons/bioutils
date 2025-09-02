# bioutils

[![Release](https://img.shields.io/github/v/release/biocommons/bioutils)](https://img.shields.io/github/v/release/biocommons/bioutils)
[![Build status](https://img.shields.io/github/actions/workflow/status/biocommons/bioutils/main.yml?branch=main)](https://github.com/biocommons/bioutils/actions/workflows/main.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/biocommons/bioutils/branch/main/graph/badge.svg)](https://codecov.io/gh/biocommons/bioutils)
[![Commit activity](https://img.shields.io/github/commit-activity/m/biocommons/bioutils)](https://img.shields.io/github/commit-activity/m/biocommons/bioutils)
[![License](https://img.shields.io/github/license/biocommons/bioutils)](https://img.shields.io/github/license/biocommons/bioutils)

Package Description

This project is a product of the [biocommons](https://biocommons.org/) community.

- **Github repository**: <https://github.com/biocommons/bioutils/>
- **Documentation** <https://biocommons.github.io/bioutils/>

## Python Package Installation

Install from PyPI with `pip install bioutils` or `uv pip install bioutils`, then try it:

## Developer Setup

### Install Prerequisites

These tools are required to get started:

- [git](https://git-scm.com/): Version control system
- [GNU make](https://www.gnu.org/software/make/): Current mechanism for consistent invocation of developer tools.
- [uv](https://docs.astral.sh/uv/): An extremely fast Python package and project manager, written in Rust.

#### MacOS or Linux Systems

- [Install brew](https://brew.sh/)
- `brew install git make uv`

#### Linux (Debian-based systems)

You may also install using distribution packages:

    sudo apt install git make

Then install uv using the [uv installation instructions](https://docs.astral.sh/uv/getting-started/installation/).

### One-time developer setup

Create a Python virtual environment, install dependencies, install pre-commit hooks, and install an editable package:

    make devready

### Development

**N.B.** Developers are strongly encouraged to use `make` to invoke tools to
ensure consistency with the CI/CD pipelines.  Type `make` to see a list of
supported targets.  A subset are listed here:

    » make
    🌟🌟 biocommons conventional make targets 🌟🌟

    Using these targets promots consistency between local development and ci/cd commands.

    usage: make [target ...]

    BASIC USAGE
    help                Display help message

    SETUP, INSTALLATION, PACKAGING
    devready            Prepare local dev env: Create virtual env, install the pre-commit hooks
    build               Build package
    publish             publish package to PyPI

    FORMATTING, TESTING, AND CODE QUALITY
    cqa                 Run code quality assessments
    test                Test the code with pytest

    DOCUMENTATION
    docs-serve          Build and serve the documentation
    docs-test           Test if documentation can be built without warnings or errors

    CLEANUP
    clean               Remove temporary and backup files
    cleaner             Remove files and directories that are easily rebuilt
    cleanest            Remove all files that can be rebuilt
    distclean           Remove untracked files and other detritus
