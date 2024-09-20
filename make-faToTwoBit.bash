#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail
set -o xtrace


userApps_version="470"
userApps_url="https://hgdownload.soe.ucsc.edu/admin/exe/userApps.archive/userApps.v${userApps_version}.src.tgz"
userApps_file="${userApps_url##*/}"
userApps_expected_checksum="f09db7e7805b1e681ca7722168eee3053ee325ad9dc8b3e0fea6b88487988260"

if [ ! -e "$userApps_file" ]; then
    curl -LOs "$userApps_url"
fi

userApps_actual_checksum=$(openssl dgst -sha256 "$userApps_file")
userApps_actual_checksum=${userApps_expected_checksum## }

if [ "$userApps_actual_checksum" != "$userApps_expected_checksum" ]; then
    echo "checksum mismatch: $userApps_file"
    exit 1
fi

tar xf "$userApps_file"
(cd userApps && patch -p1 < ../include.patch)
(cd userApps && patch -p1 < ../add-option.patch)

export PREFIX="$CONDA_PREFIX"
mkdir -p "${PREFIX}/bin"
export MACHTYPE=$(uname -m)
export BINDIR=$(pwd)/bin
export INCLUDE_PATH="${PREFIX}/include"
export LIBRARY_PATH="${PREFIX}/lib"
export LDFLAGS="${LDFLAGS:+$LDFLAGS} -L${PREFIX}/lib"
export CFLAGS="${CFLAGS:+$CFLAGS} -O3 ${LDFLAGS}"
export CXXFLAGS="${CXXFLAGS:+$CXXFLAGS} -I${PREFIX}/include ${LDFLAGS}"
export L="${LDFLAGS}"

mkdir -p "${BINDIR}"

dirs="
lib
htslib
jkOwnLib
hg/lib
utils/faToTwoBit
"

for dir in $dirs; do
    (cd "userApps/kent/src/$dir" && make)
done
cp bin/faToTwoBit "${PREFIX}/bin"
chmod 0755 "${PREFIX}/bin/faToTwoBit"

