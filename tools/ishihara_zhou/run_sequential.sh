#!/bin/bash
# Run Ishihara-Zhou cases one after another, each in its own folder (run.log
# and plot files are written there).
#
# Usage: EXE=/path/to/kynema_sgf [NP=4] [MPIRUN=mpirun] ./run_sequential.sh \
#            cases/ridge2d/smooth_klaxell_4m [more case folders ...]
# Extra ParmParse overrides can be passed through ARGS, e.g.
#   ARGS="amr.max_grid_size=64" ./run_sequential.sh ...
#
# A failed case does not stop the ones after it, but the script exits nonzero
# when any case failed, so that a caller can tell.
set -u
EXE=${EXE:?set EXE to the kynema-sgf executable}
# Each case runs from its own folder, so the executable is made an absolute
# path first: a bare name is looked up on PATH, a relative path is resolved
if [[ "$EXE" == */* ]]; then
    EXE=$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")
else
    EXE=$(command -v "$EXE") || { echo "EXE not found on PATH" >&2; exit 1; }
fi
NP=${NP:-4}
MPIRUN=${MPIRUN:-mpirun}
ARGS=${ARGS:-}
nfailed=0
for dir in "$@"; do
    echo "== $dir started $(date)"
    # shellcheck disable=SC2086
    (cd "$dir" && "$MPIRUN" -np "$NP" "$EXE" case.inp $ARGS > run.log 2>&1)
    # Kept before anything else runs, since even the $(date) below resets $?
    status=$?
    echo "== $dir finished $(date) status $status"
    if [ "$status" -ne 0 ]; then
        nfailed=$((nfailed + 1))
    fi
done
if [ "$nfailed" -ne 0 ]; then
    echo "== $nfailed of $# cases failed" >&2
    exit 1
fi
