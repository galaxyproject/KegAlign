#!/usr/bin/env bash
# Builds and runs tests/test_scoring.cpp against the real LASTZ score-set reader.
# No CUDA toolkit and no GPU are needed.
set -o errexit -o nounset -o pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# A LASTZ scoring file whose every value differs from KegAlign's built-in
# defaults, so that a value that fails to travel is visible rather than lucky.
cat > "$tmp/scores.txt" <<'SCORES'
# test score set - deliberately not HOXD70, and not the built-in bad/fill
bad_score          = X:-750
fill_score         = -55
gap_open_penalty   =   30
gap_extend_penalty =  400

     A     C     G     T
A   85  -110   -30  -120
C -110    95  -125   -30
G  -30  -125    95  -110
T -120   -30  -110    85
SCORES

cc  -std=c11   -I "$repo/common" -c -o "$tmp/scoring.o"       "$repo/common/scoring.c"
cc  -std=c11   -I "$repo/common" -c -o "$tmp/dna_utilities.o" "$repo/common/dna_utilities.c"
cc  -std=c11   -I "$repo/common" -c -o "$tmp/utilities.o"     "$repo/common/utilities.c"
c++ -std=c++17 -I "$repo/common" -c -o "$tmp/test.o"          "$repo/tests/test_scoring.cpp"
c++ -o "$tmp/test_scoring" "$tmp/test.o" "$tmp/scoring.o" "$tmp/dna_utilities.o" "$tmp/utilities.o" -lm

"$tmp/test_scoring" "$tmp/scores.txt"
