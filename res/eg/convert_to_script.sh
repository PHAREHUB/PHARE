#!/usr/bin/env bash
# usage:
#  ./res/eg/convert_to_script.sh res/eg/1d/jupyter/weak/weak.ipynb res/eg/1d/jupyter/weak/weak_notebook
#

set -ex

jupyter nbconvert --to script "$1" --output "$2"
