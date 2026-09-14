#!/bin/bash
# Run Ishihara-Zhou cases one after another, each in its own folder (run.log
# and plot files are written there).
#
# Usage: EXE=/path/to/kynema_sgf [NP=4] [MPIRUN=mpirun] ./run_sequential.sh \
#            cases/ridge2d/smooth_klaxell_4m [more case folders ...]
# Extra ParmParse overrides can be passed through ARGS, e.g.
#   ARGS="amr.max_grid_size=64" ./run_sequential.sh ...
set -u
EXE=${EXE:?set EXE to the kynema-sgf executable}
NP=${NP:-4}
MPIRUN=${MPIRUN:-mpirun}
ARGS=${ARGS:-}
for dir in "$@"; do
    echo "== $dir started $(date)"
    # shellcheck disable=SC2086
    (cd "$dir" && "$MPIRUN" -np "$NP" "$EXE" case.inp $ARGS > run.log 2>&1)
    echo "== $dir finished $(date) status $?"
done
