#!/usr/bin/env bash
# usage:
#  ./res/eg/convert_to_script.sh res/eg/1d/jupyter/weak/weak.ipynb
#

set -ex

jupyter nbconvert --to script "$1" --output run
