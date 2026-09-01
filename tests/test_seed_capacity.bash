#!/usr/bin/env bash
# Builds and runs tests/test_seed_capacity.cpp against the real common/ntcoding.cpp.
# No CUDA toolkit and no GPU are needed.
set -o errexit -o nounset -o pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

"${CXX:-c++}" -std=c++17 -Wall -I "$repo/common" \
    -o "$tmp/test_seed_capacity" \
    "$repo/tests/test_seed_capacity.cpp" \
    "$repo/common/ntcoding.cpp"

"$tmp/test_seed_capacity"
