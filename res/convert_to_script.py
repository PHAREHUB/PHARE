#!/usr/bin/env bash
# usage:
#  TODO

set -ex

jupyter nbconvert --to script my_notebook.ipynb --output my_script.py
