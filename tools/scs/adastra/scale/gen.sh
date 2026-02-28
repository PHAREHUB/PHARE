#!/usr/bin/env bash

# generates some harris scripts with increasing L0

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

set -ex

SIZ=100
XYZ=100
FOR=1
(
  cd $CWD

  mkdir -p gen

  for i in $(seq 7); do
    FILE="gen/config.$FOR.py"
    CALC=$(($XYZ*$SIZ*$SIZ *100 * 76 / 1000000000))
    echo $CALC
    echo "# L0 space requires $CALC GB" > $FILE
    cat config.top.py >> "$FILE"
    echo "cells = ($XYZ, $SIZ, $SIZ)" >> "$FILE"
    cat config.bot.py >> $FILE

    FOR=$(($FOR*2))
    XYZ=$(($XYZ*2))
  done
)

