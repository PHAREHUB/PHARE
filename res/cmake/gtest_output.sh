#!/usr/bin/env sh
# Runs a gtest binary with GTEST_OUTPUT (e.g. xml:gtest_out.xml) changed so that no two concurrent
#  processes write the same report file, which corrupts it or silently keeps only the last one:
#  - under mpirun every rank inherits GTEST_OUTPUT, only rank 0 keeps it.
#    see https://github.com/PHAREHUB/PHARE/issues/327
#  - several test binaries can share a working directory (e.g. tests/diagnostic) and ctest -j runs them
#    at once, so a report file goes in a subdirectory per binary (not named <binary>, that is the binary):
#    xml:gtest_out.xml -> xml:./gtest_reports/<binary>/gtest_out.xml
#    a report directory (xml:dir/) is kept, gtest already names the file after the binary there.
# gtest reads GTEST_OUTPUT during its own static initialization, so it is changed here, before exec.
# usage: [mpirun -n N] sh gtest_output.sh ./binary [args...]

rank=${OMPI_COMM_WORLD_RANK:-${PMIX_RANK:-${PMI_RANK:-0}}}
if [ "$rank" != "0" ]; then
  unset GTEST_OUTPUT
elif [ -n "${GTEST_OUTPUT:-}" ]; then
  binary=$(basename "$1")
  case "$GTEST_OUTPUT" in
    */) ;;
    *:*) path=${GTEST_OUTPUT#*:}
         export GTEST_OUTPUT="${GTEST_OUTPUT%%:*}:$(dirname "$path")/gtest_reports/$binary/$(basename "$path")" ;;
    *) export GTEST_OUTPUT="$GTEST_OUTPUT:gtest_reports/$binary/" ;;
  esac
fi
exec "$@"
