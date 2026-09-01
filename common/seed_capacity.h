#pragma once

#include <stdint.h>

// The largest number of seeds a single call to SeedAndFilter() can be handed,
// and therefore the size that d_seed_offsets (and d_hit_num_vec) must be
// allocated to on the GPU.
//
// src/seeder.cpp walks at most wga_chunk_size query positions per call, and for
// each position whose k-mer is valid it pushes one seed for the position itself
// plus one more for every position in the seed pattern at which a transition is
// allowed.  This function has to stay in step with that loop; tests/ pins the
// two together.
static inline uint64_t MaxSeedsPerChunk(bool transition, uint32_t num_transitions, uint32_t wga_chunk_size)
{
    if (!transition) {
        return (uint64_t) wga_chunk_size;
    }

    return (1ULL + num_transitions) * wga_chunk_size;
}
