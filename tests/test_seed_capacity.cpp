// Tests that the GPU seed buffer is sized for the seed pattern actually in use.
//
// Run with:  bash tests/test_seed_capacity.bash
//
// Needs a C++ compiler, but no CUDA toolkit and no GPU: the rule under test is
// MaxSeedsPerChunk(), which is pure arithmetic, and the transition counting is
// done by the real common/ntcoding.cpp.
//
// ---------------------------------------------------------------------------
// The bug this all exists for: MAX_SEEDS was a hardcoded 13*wga_chunk_size,
// where 13 is 1 + the 12 transition positions of the DEFAULT 12of19 seed. The
// 14of22 seed has 14 transition positions and so emits 15 seeds per query
// position, overrunning d_seed_offsets. On a real usegalaxy.org run this
// surfaced as
//
//     Error: cudaMemcpy of 26650680 bytes for seed_offsets failed
//            with error " invalid argument "
//
// 26650680/8 = 3331335 seeds copied into a buffer sized for 3250000. The
// assert() that should have caught it is compiled out by NDEBUG in the Release
// builds conda-forge and Galaxy ship.
// ---------------------------------------------------------------------------

#include <stdio.h>
#include <stdint.h>
#include <string>

#include "ntcoding.h"
#include "seed_capacity.h"

static int failures = 0;

static void expect_eq(const char* label, uint64_t actual, uint64_t expected)
{
    if (actual != expected) {
        printf("not ok - %s\n     got: %lu\nexpected: %lu\n", label, actual, expected);
        failures++;
    } else {
        printf("ok - %s\n", label);
    }
}

static void expect_at_least(const char* label, uint64_t actual, uint64_t floor)
{
    if (actual < floor) {
        printf("not ok - %s\n     got: %lu\nexpected at least: %lu (short by %lu)\n",
               label, actual, floor, floor - actual);
        failures++;
    } else {
        printf("ok - %s\n", label);
    }
}

// Mirrors the push sequence in src/seeder.cpp: for a query position whose k-mer
// is valid, one seed for the position itself, then one more for each position in
// the pattern at which a transition is allowed.  Uses the real IsTransitionAtPos()
// so that a change to the seed parser shows up here.
static uint64_t SeedsPushedPerPosition(int kmer_size)
{
    uint64_t pushed = 1;
    for (int t = 0; t < kmer_size; t++) {
        if (IsTransitionAtPos(t) == 1) {
            pushed++;
        }
    }
    return pushed;
}

// The seeder hands SeedAndFilter() at most wga_chunk_size positions at a time.
static uint64_t WorstCaseSeedsInOneChunk(const std::string& shape, uint32_t wga_chunk_size)
{
    int kmer_size = GenerateShapePos(shape);
    return SeedsPushedPerPosition(kmer_size) * wga_chunk_size;
}

// main.cpp normalises a user's pattern so every match position becomes 'T'
// before calling GenerateShapePos(), so these are the shapes as the parser sees them.
static const char* SEED_12of19 = "TTT0T00TT00T0T0TTTT";
static const char* SEED_14of22 = "TTT0T0TT00TT00T0T0TTTT";

int main()
{
    const uint32_t chunk = 250000;   // DEFAULT_WGA_CHUNK, and the failing job's value

    // --- the default seed must be unaffected -------------------------------
    GenerateShapePos(SEED_12of19);
    expect_eq("12of19 has 12 transition positions", GetNumTransitions(), 12);
    expect_eq("12of19 capacity is unchanged at 13*chunk",
              MaxSeedsPerChunk(true, GetNumTransitions(), chunk), 13ULL * chunk);
    expect_at_least("12of19 capacity covers what the seeder pushes",
                    MaxSeedsPerChunk(true, GetNumTransitions(), chunk),
                    WorstCaseSeedsInOneChunk(SEED_12of19, chunk));

    // --- the bug -----------------------------------------------------------
    GenerateShapePos(SEED_14of22);
    expect_eq("14of22 has 14 transition positions", GetNumTransitions(), 14);
    expect_at_least("14of22 capacity covers what the seeder pushes",
                    MaxSeedsPerChunk(true, GetNumTransitions(), chunk),
                    WorstCaseSeedsInOneChunk(SEED_14of22, chunk));

    // The seed count from the job that actually died.
    GenerateShapePos(SEED_14of22);
    expect_at_least("14of22 capacity covers the 3331335 seeds of the failing run",
                    MaxSeedsPerChunk(true, GetNumTransitions(), chunk), 3331335);

    // --- sized by transitions, not by weight -------------------------------
    // Mixed pattern: '1' is a match position that does NOT allow a transition,
    // so weight is 12 but only 11 transitions -> 12 seeds per position. A fix
    // that reached for kmer_size instead of the transition count would over-allocate.
    const char* mixed = "TTT0T00TT00T0T0TT1T";
    GenerateShapePos(mixed);
    expect_eq("mixed pattern: weight 12 but only 11 transitions", GetNumTransitions(), 11);
    expect_eq("mixed pattern capacity follows transitions, not weight",
              MaxSeedsPerChunk(true, GetNumTransitions(), chunk),
              WorstCaseSeedsInOneChunk(mixed, chunk));

    // --- the count must not carry over between patterns ---------------------
    GenerateShapePos(SEED_14of22);
    GenerateShapePos(SEED_12of19);
    expect_eq("GenerateShapePos resets the transition count", GetNumTransitions(), 12);

    // --- --notransition path ------------------------------------------------
    GenerateShapePos(SEED_14of22);
    expect_eq("notransition gives one seed per position",
              MaxSeedsPerChunk(false, GetNumTransitions(), chunk), chunk);

    if (failures > 0) {
        printf("\n%d test(s) failed\n", failures);
        return 1;
    }
    printf("\nall tests passed\n");
    return 0;
}
