#!/usr/bin/env sh
# Runs a gtest binary under mpirun with only rank 0 writing GTEST_OUTPUT (e.g. xml:gtest_out.xml),
#  every rank inherits it and writing the same file concurrently corrupts it.
#  see https://github.com/PHAREHUB/PHARE/issues/327
# gtest reads GTEST_OUTPUT during its own static initialization, so it is unset here, before exec.
# usage: mpirun -n N sh gtest_rank0_output.sh ./binary [args...]

rank=${OMPI_COMM_WORLD_RANK:-${PMIX_RANK:-${PMI_RANK:-0}}}
[ "$rank" = "0" ] || unset GTEST_OUTPUT
exec "$@"
