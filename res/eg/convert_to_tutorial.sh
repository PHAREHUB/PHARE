#!/usr/bin/env bash
# usage:
#  ./res/eg/convert_to_tutorial.sh res/eg/1d/jupyter/weak/weak.ipynb
#

set -ex

jupyter nbconvert --to notebook "$1"  \
   --TagRemovePreprocessor.enabled=True            \
   --TagRemovePreprocessor.remove_cell_tags answer  \
   --output tutorial
