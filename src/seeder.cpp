#include <algorithm>
#include <numeric>
#include <tuple>
#include <vector>

#include "graph.h"
#include "ntcoding.h"
#include "seed_filter.h"
#include "store.h"

std::atomic<uint64_t> seeder_body::num_seed_hits(0);
std::atomic<uint64_t> seeder_body::num_seeds(0);
std::atomic<uint64_t> seeder_body::num_hsps(0);
std::atomic<uint32_t> seeder_body::total_xdrop(0);
std::atomic<uint32_t> seeder_body::num_seeded_regions[BUFFER_DEPTH]={};


// Remove duplicate HSPs once the chunks of this interval have been joined.
//
// SeedAndFilter() does sort and unique_copy() its output, but once per
// MAX_HITS-bounded iteration, concatenating the results -- so duplicates already
// survive within a single chunk. seeder.cpp then concatenates chunks with no
// dedup at all. An HSP reachable from seeds in two places is emitted twice, and
// lastz repeats the same gapped extension and reports a duplicate alignment.
//
// Returns how many were removed, so the #HSPs counter stays consistent with what
// actually reaches the segments files.
//
// stable_sort, not sort: the survivor of a group has to be the FIRST occurrence.
// With an unstable sort, which copy is kept -- and therefore the byte content of
// the emitted file -- depends on the standard library's unspecified ordering of
// equal elements, and is not reproducible across toolchains.
//
// Compaction is in place. Building a filtered copy would roughly double peak
// memory for this call, and seeder runs tbb::flow::unlimited across every core.
//
// SCOPE, stated plainly. This removes bit-identical HSPs within one interval.
// It deliberately does NOT reproduce the GPU's hspEqual(), which also collapses
// an HSP contained by another on the same diagonal: applying containment here
// could drop an HSP that lastz would have reported. It also does not reach
// across intervals, which are separate calls writing separate files.
static size_t drop_duplicate_hsps (std::vector<segmentPair>& hsps) {

    if (hsps.size() < 2) {
        return 0;
    }

    auto key = [&hsps](uint32_t i) {
        const segmentPair& s = hsps[i];
        return std::make_tuple(s.ref_start, s.query_start, s.len, s.score);
    };

    std::vector<uint32_t> order(hsps.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(),
                     [&key](uint32_t a, uint32_t b) { return key(a) < key(b); });

    std::vector<bool> duplicate(hsps.size(), false);
    size_t removed = 0;
    for (size_t i = 1; i < order.size(); i++) {
        if (key(order[i]) == key(order[i-1])) {
            duplicate[order[i]] = true;
            removed++;
        }
    }

    if (removed == 0) {
        return 0;
    }

    size_t keep = 0;
    for (size_t i = 0; i < hsps.size(); i++) {
        if (!duplicate[i]) {
            hsps[keep++] = hsps[i];
        }
    }
    hsps.resize(keep);

    return removed;
}

printer_input seeder_body::operator()(seeder_input input) {

    auto &payload = std::get<0>(input);
    size_t token  = std::get<1>(input);

    auto &block_data    = std::get<0>(payload);
    auto &interval_data = std::get<1>(payload);

    int r_block_index    = block_data.r_index;
    int q_block_index    = block_data.q_index;
    size_t r_block_start = block_data.r_start;
    size_t q_block_start = block_data.q_start;
    uint32_t r_block_len = block_data.r_len;
    uint32_t q_block_len = block_data.q_len;

    uint32_t q_inter_start = interval_data.start;
    uint32_t q_inter_end   = interval_data.end;
    uint32_t buffer        = interval_data.buffer;
    uint32_t num_invoked   = interval_data.num_invoked;
    uint32_t num_intervals = interval_data.num_intervals;

    uint32_t rc_q_inter_start = q_block_len - q_inter_end;
    uint32_t rc_q_inter_end   = q_block_len - q_inter_start;

    uint64_t kmer_index;
    uint64_t transition_index;
    uint64_t seed_offset;

    std::vector<segmentPair> fw_hsps;
    std::vector<segmentPair> rc_hsps;
    fw_hsps.clear();
    rc_hsps.clear();

    fprintf (stderr, "Query block %u, interval %u/%u (%u:%u) with buffer %u\n", q_block_index, num_invoked, num_intervals, q_inter_start, q_inter_end, buffer);

    if(cfg.strand == "plus" || cfg.strand == "both"){
        for (uint32_t i = q_inter_start; i < q_inter_end; i += cfg.wga_chunk_size) {

            //end position
            uint32_t e = std::min(i + cfg.wga_chunk_size, q_inter_end + 1);

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();

            //start to end position in the chunk
            for (uint32_t j = i; j < e; j++) {

                kmer_index = GetKmerIndexAtPos(query_DRAM->buffer, q_block_start+j, cfg.seed.size);
                if (kmer_index != ((uint32_t) 1 << 31)) {
                    seed_offset = (kmer_index << 32) + j;
                    seed_offset_vector.push_back(seed_offset); 

                    if (cfg.seed.transition) {
                        for (int t=0; t < cfg.seed.kmer_size; t++) {
                            if (IsTransitionAtPos(t) == 1) {
                                transition_index = (kmer_index ^ (TRANSITION_MASK << (2*t)));
                                seed_offset = (transition_index << 32) + j;
                                seed_offset_vector.push_back(seed_offset); 
                            }
                        }
                    }
                }
            }

            if(seed_offset_vector.size() > 0){
                seeder_body::num_seeds += seed_offset_vector.size();
                std::vector<segmentPair> anchors = g_SeedAndFilter(seed_offset_vector, false, buffer);
                seeder_body::num_seed_hits += anchors[0].score;
                if(anchors.size() > 1){
                    fw_hsps.insert(fw_hsps.end(), anchors.begin()+1, anchors.end());
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }
    }

    if(cfg.strand == "minus" || cfg.strand == "both"){
        for (uint32_t i = rc_q_inter_start; i < rc_q_inter_end; i += cfg.wga_chunk_size) {
            uint32_t e = std::min(i + cfg.wga_chunk_size, rc_q_inter_end + 1);

            std::vector<uint64_t> seed_offset_vector;
            seed_offset_vector.clear();
            for (uint32_t j = i; j < e; j++) {
                kmer_index = GetKmerIndexAtPos(query_rc_DRAM->buffer, q_block_start+j, cfg.seed.size);
                if (kmer_index != ((uint32_t) 1 << 31)) {
                    seed_offset = (kmer_index << 32) + j;
                    seed_offset_vector.push_back(seed_offset); 
                    if (cfg.seed.transition) {
                        for (int t=0; t < cfg.seed.kmer_size; t++) {
                            if (IsTransitionAtPos(t) == 1) {
                                transition_index = (kmer_index ^ (TRANSITION_MASK << (2*t)));
                                seed_offset = (transition_index << 32) + j;
                                seed_offset_vector.push_back(seed_offset); 
                            }
                        }
                    }
                }
            }

            if(seed_offset_vector.size() > 0){
                seeder_body::num_seeds += seed_offset_vector.size();
                std::vector<segmentPair> anchors = g_SeedAndFilter(seed_offset_vector, true, buffer);
                seeder_body::num_seed_hits += anchors[0].score;
                if(anchors.size() > 1){
                    rc_hsps.insert(rc_hsps.end(), anchors.begin()+1, anchors.end());
                    seeder_body::num_hsps += anchors.size()-1;
                }
            }
        }
    }

    // Adjust the counter too: it is incremented per chunk above, so without this
    // the #HSPs figure --debug prints would exceed what the segments files hold.
    seeder_body::num_hsps -= drop_duplicate_hsps(fw_hsps);
    seeder_body::num_hsps -= drop_duplicate_hsps(rc_hsps);

    seeder_body::num_seeded_regions[buffer] += 1;
    seeder_body::total_xdrop += 1;

    return printer_input(printer_payload(payload, fw_hsps, rc_hsps), token);
}
