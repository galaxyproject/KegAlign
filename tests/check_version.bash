#!/usr/bin/env bash
# Checks that the VERSION define in common/parameters.h matches an expected
# version string -- in CI, the tag being released.
#
# Run with:  bash tests/check_version.bash v0.2.2.14
#
# This exists because #define VERSION sat at v0.1.2.8 from that release all the
# way through v0.2.1.13: five releases reported the wrong version from
# --version and in the usage banner, because nothing tied the define to the tag.
set -o errexit -o nounset -o pipefail

if [ "$#" -ne 1 ]; then
    echo "usage: ${0##*/} <expected-version>" >&2
    exit 2
fi

expected="$1"
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
header="$repo/common/parameters.h"

defined="$(sed -n 's/^#define VERSION "\(.*\)"$/\1/p' "$header")"

if [ -z "$defined" ]; then
    echo "not ok - no VERSION define found in $header" >&2
    exit 1
fi

if [ "$defined" != "$expected" ]; then
    echo "not ok - $header defines VERSION $defined, but the tag is $expected" >&2
    exit 1
fi

echo "ok - VERSION is $defined"
