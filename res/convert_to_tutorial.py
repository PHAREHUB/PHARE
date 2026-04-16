#!/usr/bin/env bash
# usage:
#  TODO

jupyter nbconvert --to notebook my_notebook.ipynb  \
   --TagRemovePreprocessor.enabled=True            \
   -TagRemovePreprocessor.remove_cell_tags answer  \
   --output my_notebook_without_answer.ipynb
