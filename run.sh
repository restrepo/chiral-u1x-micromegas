#!/usr/bin/env bash
# Run micrOMEGAs for the U1X model on a parameter file.
# Usage:  ./run.sh benchmarks/BP1.par
set -euo pipefail

REPO=$(cd "$(dirname "$0")" && pwd)
BUILD=${BUILD:-$REPO/build}
PROJ=$BUILD/micromegas_7.1.1/U1X

[[ $# -eq 1 ]] || { echo "Usage: $0 <file.par>"; exit 1; }
[[ -x "$PROJ/main" ]] || { echo "Run ./setup.sh first"; exit 1; }

PAR=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
cd "$PROJ"
./main "$PAR" 2>&1 | grep -v '^PROCESS\|Processes of this type are absent\|^Delete diagrams'
