#!/usr/bin/env bash
#
# PHLOP_SCOPE timer comparison of two PHARE builds, eg master vs a PR
#
#   tools/bench/functional/compare.sh <base checkout> <new checkout> [scope_compare args]
#   tools/bench/functional/compare.sh ../master . -t instructions=1 -t time=5
#
# Each checkout runs its own pyphare and build (env BUILD, default build, relative
# to the checkout), both run this script's bench, so the config is the same.
# Runs go to $OUT/{base,new}/runs/<config>/repeat=<i>/, then phlop.timing.scope_compare
# compares each config, exit 1 if any regressed. Each repeat is its own process, base
# and new alternating, so both start from the same state and see the same machine.
#
#  env OUT: default .phare_bench/compare
#  env NP: mpirun -n NP, default serial. Both sides run the same NP
#  env MPIRUN: the launcher, default "mpirun -n $NP", eg "mpirun -n 4 --bind-to core"
#  env LAUNCH: prefix of python3, after MPIRUN, eg "setarch -R" for no ASLR
#  env BENCH: default compare_harris_2d.py next to this script
#  env REPEATS: runs per side, default 3
#  env PHARE_BENCH_STEPS: see compare_harris_2d.py
#  env PHLOP_PERF_ALL / PHLOP_PERF_LIGHT: perf counters, see phlop
#  env PHLOP_SCOPE_MEM=1 (or peak): memory of mem scopes, gate with eg -t mem_peak=2

set -eu

CWD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
(( $# >= 2 )) || { sed -n '3,22p' "${BASH_SOURCE[0]}"; exit 2; }

BASE="$(cd "$1" && pwd)"
NEW="$(cd "$2" && pwd)"
shift 2
BUILD="${BUILD:-build}"
OUT="$(mkdir -p "${OUT:-.phare_bench/compare}" && cd "${OUT:-.phare_bench/compare}" && pwd)"
BENCH="${BENCH:-$CWD/compare_harris_2d.py}"
REPEATS="${REPEATS:-3}"
RUN=() PREFIX=()
[ -n "${NP:-}" ] && RUN=(mpirun -n "$NP")
[ -n "${MPIRUN:-}" ] && read -ra RUN <<< "$MPIRUN"
[ -n "${LAUNCH:-}" ] && read -ra PREFIX <<< "$LAUNCH"

run() { # <label> <checkout> <repeat>
  local out="$OUT/$1"
  (
    cd "$out/work" # each side its own .phare and .log
    PYTHONPATH="$2:$2/$BUILD:$2/pyphare${PYTHONPATH:+:$PYTHONPATH}" \
      PHARE_SCOPE_TIMING=1 PHARE_BENCH_TAR=0 PHARE_BENCH_MON=0 PHARE_BENCH_OUT="$out/runs" \
      PHARE_BENCH_REPEAT="$3" "${RUN[@]}" "${PREFIX[@]}" python3 -Ou "$BENCH"
  )
}

for side in base new; do
  rm -rf "${OUT:?}/$side"
  mkdir -p "$OUT/$side/work"
done
for (( i = 0; i < REPEATS; i++ )); do
  echo "== repeat $i"
  run base "$BASE" "$i"
  run new "$NEW" "$i"
done

ec=0
for cfg in "$OUT"/base/runs/*/; do
  cfg="$(basename "$cfg")"
  echo "== $cfg"
  if [ ! -d "$OUT/new/runs/$cfg" ]; then
    echo "missing in new" && ec=2 && continue
  fi
  e=0
  python3 -m phlop.timing.scope_compare "$OUT/base/runs/$cfg" "$OUT/new/runs/$cfg" \
    -p 'rank.*.txt' "$@" || e=$?
  if (( e > ec )); then ec=$e; fi
done
exit $ec
