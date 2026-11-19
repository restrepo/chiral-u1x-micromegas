#!/usr/bin/env bash
# Download and build micrOMEGAs 7.1.1 and LanHEP 4.0.0, generate the U1X model
# with LanHEP and compile the micrOMEGAs executable.
#
# Usage:  ./setup.sh [--no-lanhep]
#   --no-lanhep  use the CalcHEP model files shipped in model/calchep instead of
#                running LanHEP (LanHEP is then not downloaded).
# Everything is built under ./build (override with BUILD=/some/dir ./setup.sh).

set -euo pipefail

REPO=$(cd "$(dirname "$0")" && pwd)
BUILD=${BUILD:-$REPO/build}

MO_VERSION=7.1.1
MO_TGZ=micromegas_$MO_VERSION.tgz
MO_URL=https://zenodo.org/api/records/20547050/files/$MO_TGZ/content
MO_MD5=291ce2a028c408a2b1516655e18359ea

LH_TGZ=lhep400.tgz
LH_URL=https://theory.sinp.msu.ru/~semenov/$LH_TGZ
LH_MD5=239d5675f57b08192c04199ced133057

USE_LANHEP=1
[[ "${1:-}" == "--no-lanhep" ]] && USE_LANHEP=0

download() {  # url file md5
    if [[ ! -f "$2" ]]; then
        echo ">>> downloading $2"
        curl -fSL -o "$2" "$1"
    fi
    echo "$3  $2" | md5sum -c -
}

mkdir -p "$BUILD/downloads"
cd "$BUILD/downloads"
download "$MO_URL" "$MO_TGZ" "$MO_MD5"
[[ $USE_LANHEP == 1 ]] && download "$LH_URL" "$LH_TGZ" "$LH_MD5"

cd "$BUILD"
MO=$BUILD/micromegas_$MO_VERSION
if [[ ! -f "$MO/lib/micromegas.a" ]]; then
    [[ -d "$MO" ]] || tar xzf "downloads/$MO_TGZ"
    echo ">>> building micrOMEGAs (a few minutes, log in $MO/make.log)"
    make -C "$MO" > "$MO/make.log" 2>&1
fi

if [[ $USE_LANHEP == 1 ]]; then
    LH=$BUILD/lanhep400
    if [[ ! -x "$LH/lhep" ]]; then
        [[ -d "$LH" ]] || tar xzf "downloads/$LH_TGZ"
        echo ">>> building LanHEP (log in $LH/make.log)"
        make -C "$LH" > "$LH/make.log" 2>&1
    fi
fi

PROJ=$MO/U1X
if [[ ! -d "$PROJ" ]]; then
    echo ">>> creating micrOMEGAs project U1X"
    (cd "$MO" && ./newProject U1X > /dev/null)
fi

mkdir -p "$PROJ/work/lanhep"
cp "$REPO"/model/lanhep/*.lhp "$REPO"/model/lanhep/*.inc "$PROJ/work/lanhep/"
if [[ $USE_LANHEP == 1 ]]; then
    echo ">>> running LanHEP"
    (cd "$PROJ/work/lanhep" && "$LH/lhep" -ca source.lhp | grep -v processed)
    cp "$PROJ"/work/lanhep/{prtcls1,vars1,func1,lgrng1}.mdl "$PROJ/work/models/"
else
    cp "$REPO"/model/calchep/{prtcls1,vars1,func1,lgrng1}.mdl "$PROJ/work/models/"
fi

cp "$REPO/main.c" "$PROJ/main.c"
echo ">>> compiling main.c"
make -C "$PROJ" main=main.c > "$PROJ/build.log" 2>&1

echo ">>> done. Run a benchmark with:  ./run.sh benchmarks/BP1.par"
