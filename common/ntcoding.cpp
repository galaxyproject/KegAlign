#include <assert.h>
#include <string>
#include "ntcoding.h"
#include "parameters.h"

int shape_pos[32];
int shape_size;
int transition_pos[32];
int num_transitions;

inline uint32_t NtChar2Int (char nt) {
    switch(nt) {
        case 'A': return A_NT;
        case 'C': return C_NT;
        case 'G': return G_NT;
        case 'T': return T_NT;
        case 'N': return N_NT;
        default: return N_NT;
    }
}

// Returns the pattern's weight -- its true count of match positions, which may
// exceed what shape_pos[] can hold.  The caller is expected to reject anything
// above MAX_SEED_WEIGHT; this function must not corrupt memory while producing
// the number that check is made against, so it stops writing at the bound and
// keeps counting.
int GenerateShapePos (std::string shape) {
    shape_size = 0;
    num_transitions = 0;
    int j = 0;
    int weight = 0;
    for (size_t i = 0; i < shape.length(); i++) {
        if ((shape[i] == '1') || (shape[i] == 'T')) {
            weight++;
            if (shape_size >= MAX_SEED_WEIGHT) {
                continue;
            }
            shape_pos[shape_size++] = i;
            if (shape[i] == 'T') {
                transition_pos[j] = 1;
                num_transitions++;
            }
            else {
                transition_pos[j] = 0;
            }
            j++;
        }
    }
    return weight;
}

int IsTransitionAtPos(int t) {
    return transition_pos[t];
}

// How many positions in the current seed pattern allow a transition.  This is
// what bounds the size of a seed_offset_vector, so MaxSeedsPerChunk() needs it.
int GetNumTransitions() {
    return num_transitions;
}

uint32_t GetKmerIndexAtPos (char* sequence, size_t pos, uint32_t seed_size) {

    uint32_t nt[MAX_SEED_SPAN];

    for(uint32_t i = 0; i < seed_size; i++){
        nt[i] = NtChar2Int(sequence[pos+i]);
        if (nt[i] == N_NT) {
            return INVALID_KMER;
        }
    }

    uint32_t kmer = 0;

    for (int i = 0; i < shape_size; i++) {
        kmer = (kmer << 2) + nt[shape_pos[i]];
    }

    return kmer;
}

const char *rev_comp =
    "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
    "NNNNNN&NNNNNNNNNNNNNNNNNNNNNNNNN"
    "NTVGHNNCDNNMNKNNNNYSAABWNRNNNNNN"
    "NtvghnncdnnmnknnnnysaabwnrnNNNNN";

void RevComp(char* dst_buffer, char* src_buffer, size_t rc_start, size_t start, size_t len){

    size_t r = rc_start;
    for (size_t i = start+len; i> start; i--) {
        int rev_comp_index = static_cast<int>(src_buffer[i-1]);
        if (rev_comp_index >= 0 && rev_comp_index <= 127) {
            dst_buffer[r++] = rev_comp[rev_comp_index];
        } else {
            dst_buffer[r++] = 'N';
            fprintf(stderr, "Bad Nt char! '%c' %lu\n", src_buffer[i-1], i-1);
        }
    }
}
