#!/usr/bin/env bash
# Differential parity test: KegAlign's gap-free extension against LASTZ's.
#
# Run with:  bash tests/test_lastz_parity.bash [target.fa query.fa]
#
# KegAlign reimplements LASTZ's seed-filter-extend paradigm on the GPU, then hands
# the resulting HSPs to LASTZ (--segments=) for gapped extension. Every other test
# in tests/ checks a piece of that reimplementation against its own arithmetic.
# This one checks it against the original: run both tools over the same pair with
# the same parameters, ask each for its HSPs in the same eight columns, and diff.
#
# It runs TWO cases. The second soft-masks both inputs, because masked bases are
# the interesting ones: they are excluded from seeding and score bad_score in
# extension, and they are the only way to reach the L_NT/N_NT/X_NT/E_NT paths in
# find_hsps. Measured on the default pair: masking raises the share of HSPs inside
# the [hspthresh, 3*hspthresh] entropy band from 35% to 55%.
#
# Masks are applied at identical FORWARD coordinates, so target and query masked
# regions correspond on the plus strand only; the query is reverse-complemented
# for the minus strand. The masked case still exercises both strands, but the
# aligned-masked-pair rationale holds for + only.
#
# NOTE what this does NOT reach: masked columns score bad_score against the
# default xdrop, so they x-drop immediately and never land inside a counted
# extent. Exercising the find_hsps entropy counters needs IUPAC codes in the
# input, or --ambiguous=n, not soft-masking.
#
# ⚠ NEEDS A GPU, plus `kegalign` and `lastz` on PATH. Everything else in tests/ is
# deliberately host-side and GPU-free so CI can run it; this one cannot be, so it
# is a manual pre-release check. Exits 77 (skip) when it cannot run.
#
# WHY IT EXISTS. Reading half of a scoring-file reader is not evidence about
# behaviour. On 2026-09-02 a "lastz masking parity" feature was written, compiled
# and tested on the strength of one such reading, and was wrong -- LASTZ builds two
# score sets (see masked_score_set() in common/dna_utilities.c) and the ungapped
# stage uses the masked one, so KegAlign was already at parity. A five-minute run
# of this test would have refuted it. Prefer this over source archaeology.
#
# FIELD NAMES ARE NOT GUESSABLE. The general-format column is "strand2", not
# "strand" -- LASTZ rejects the latter. Run `lastz ... --format=general` with no
# field list to have it print its own header, which is authoritative. Note zstart*
# is 0-based half-open while start* is 1-based closed; KegAlign writes 1-based.
#
# The defaults are the repo's own test-data: 31 sequences, ~372 kb and ~377 kb,
# which fit in ONE sequence block (default seq_block_size is 500 Mb). That matters.
# KegAlign separates blocks with '&' and will not extend an HSP across one, while
# LASTZ has no such boundary, so on a multi-block input some disagreement near
# every block edge is EXPECTED and is not a parity defect.
set -o errexit -o nounset -o pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# Held identical across both tools: KegAlign's defaults, spelled out.
SEED=12of19 ; XDROP=910 ; HSPTHRESH=3000 ; STEP=1

need() { command -v "$1" >/dev/null || { echo "SKIP: $1 not on PATH"; exit 77; }; }
need kegalign
need lastz
command -v nvidia-smi >/dev/null && nvidia-smi -L >/dev/null 2>&1 \
    || { echo "SKIP: no usable GPU (kegalign cannot run)"; exit 77; }

if [ "$#" -eq 2 ]; then
    cp "$1" "$tmp/target.fa" ; cp "$2" "$tmp/query.fa"
else
    gzip -cdf "$repo/test-data/apple.fasta.gz"  > "$tmp/target.fa"
    gzip -cdf "$repo/test-data/orange.fasta.gz" > "$tmp/query.fa"
fi

# Soft-mask a deterministic 50%: lower-case the second half of every 1000 bp
# window. Applied at IDENTICAL coordinates in both files so that homologous masked
# regions align, which is what produces aligned masked pairs during extension.
soft_mask() {
    # Lower-case the second half of every 1000 bp window, counted per sequence.
    # Emits run-at-a-time with printf rather than accumulating a string: "out = out c"
    # per character is quadratic in line length and takes ~18 minutes on a 10 Mb
    # single-line record, which one-record-per-line FASTA routinely has.
    # Existing lower case is left alone rather than upper-cased, so a repeat-masked
    # input is not silently unmasked -- which does mean the masked fraction of such
    # an input exceeds 50%.
    awk '
        /^>/ { print; pos = 0; next }
        {
            n = length($0); i = 1
            while (i <= n) {
                p   = pos % 1000
                run = (p < 500) ? 500 - p : 1000 - p
                if (run > n - i + 1) run = n - i + 1
                seg = substr($0, i, run)
                printf "%s", (p < 500) ? seg : tolower(seg)
                i += run; pos += run
            }
            printf "\n"
        }' "$1" > "$2"
}
soft_mask "$tmp/target.fa" "$tmp/target.masked.fa"
soft_mask "$tmp/query.fa"  "$tmp/query.masked.fa"

failures=0

run_case() {
    local label="$1" target="$2" query="$3" dir="$tmp/$1"
    mkdir -p "$dir/work"

    # Guarded: under errexit a bare lastz failure would tear down the whole script
    # mid-case instead of recording one failure and running the other case.
    lastz "$target[multiple]" "$query" \
        --nogapped --strand=both \
        --seed="$SEED" --xdrop="$XDROP" --hspthresh="$HSPTHRESH" --step="$STEP" \
        --format=general:name1,start1,end1,name2,start2,end2,strand2,score \
        --output="$dir/lastz.raw" \
        || { echo "not ok - $label: lastz failed"; failures=$((failures+1)); return; }
    grep -v '^#' "$dir/lastz.raw" | sort > "$dir/lastz.hsp"

    ( cd "$dir" && kegalign "$target" "$query" work/ \
        --seed "$SEED" --xdrop "$XDROP" --hspthresh "$HSPTHRESH" --step "$STEP" \
        --strand both > commands.txt 2> kegalign.err ) \
        || { echo "not ok - $label: kegalign failed"; tail -3 "$dir/kegalign.err"; failures=$((failures+1)); return; }
    cat "$dir"/*.segments 2>/dev/null | sort > "$dir/kegalign.hsp" || : > "$dir/kegalign.hsp"

    local l k both dup
    l=$(wc -l < "$dir/lastz.hsp") ; k=$(wc -l < "$dir/kegalign.hsp")
    both=$(comm -12 "$dir/lastz.hsp" "$dir/kegalign.hsp" | wc -l)
    # Per FILE, not across the concatenation: dedup happens per interval, and each
    # interval writes its own .segments consumed by its own lastz call. Counting
    # across files would flag cross-interval duplicates, which the fix does not
    # claim to remove, on any query longer than lastz_interval_size.
    dup=0
    for f in "$dir"/*.segments; do
        [ -e "$f" ] || continue
        # total - distinct, not `uniq -d | wc -l`: that counts duplicated GROUPS,
        # so a segment emitted three times would report as one.
        dup=$(( dup + $(wc -l < "$f") - $(sort -u "$f" | wc -l) ))
    done
    printf "  %-8s lastz=%-6s kegalign=%-6s identical=%-6s duplicates=%s\n" \
           "$label" "$l" "$k" "$both" "$dup"

    # Both directions. "both == l" alone only asserts KegAlign is a SUPERSET, so a
    # regression emitting spurious distinct HSPs would pass silently.
    if [ "$both" -ne "$l" ] || [ "$k" -ne "$l" ]; then
        echo "not ok - $label: HSP sets differ (missing $((l-both)), surplus $((k-both)))"
        # A systematic coordinate-convention mismatch makes EVERY line differ, which
        # reads like a catastrophe rather than an off-by-one. Rule that out first.
        local s
        s=$(awk 'BEGIN{FS=OFS="\t"} {$3=$3-1; $6=$6-1; print}' "$dir/kegalign.hsp" \
            | sort | comm -12 "$dir/lastz.hsp" - | wc -l)
        [ "$s" -gt "$both" ] && echo "    NOTE: shifting KegAlign end coords by -1 matches $s (vs $both): systematic off-by-one"
        echo "    LASTZ-only:";    comm -23 "$dir/lastz.hsp" "$dir/kegalign.hsp" | head -3 | sed 's/^/      /'
        echo "    KegAlign-only:"; comm -13 "$dir/lastz.hsp" "$dir/kegalign.hsp" | head -3 | sed 's/^/      /'
        failures=$((failures+1))
    elif [ "$dup" -ne 0 ]; then
        echo "not ok - $label: $dup duplicate segment(s); lastz would repeat that gapped extension"
        sort "$dir/kegalign.hsp" | uniq -d | head -3 | sed 's/^/      /'
        failures=$((failures+1))
    else
        echo "ok - $label: KegAlign and LASTZ agree exactly, no duplicates"
    fi
}

run_case unmasked "$tmp/target.fa"        "$tmp/query.fa"
run_case masked   "$tmp/target.masked.fa" "$tmp/query.masked.fa"

if [ "$failures" -gt 0 ]; then
    echo ; echo "$failures case(s) failed"
    exit 1
fi
echo ; echo "all cases passed"
