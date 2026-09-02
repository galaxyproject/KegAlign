#define VERSION "v0.3.0"

// Upper bound on a seed pattern's weight (its count of match positions).
//
// 15, not 32, and the tighter bound is the load-bearing one: GetKmerIndexAtPos
// packs two bits per match position into a uint32_t, and INVALID_KMER is 1<<31.
// At weight 16 a legitimate k-mer can equal that sentinel and be discarded as
// invalid; above 16 the k-mer wraps outright.  Staying at or below 15 also keeps
// shape_pos[32]/transition_pos[32] in bounds, and keeps GetNumTransitions()
// honest -- MaxSeedsPerChunk() sizes the GPU seed buffer from that count, so a
// corrupted one silently breaks the capacity invariant.
#define MAX_SEED_WEIGHT 15

// Lower bound on a seed pattern's weight.
//
// This is the bound GenerateSeedPosTable() already states -- assert(kmer_size > 3)
// in common/seed_pos_table.cu -- promoted to something the shipped binary
// enforces. NDEBUG deletes that assert in a Release build, the same way it
// deleted the seed-buffer assert, which is how the 14of22 overflow reached a
// user as a bare cudaMemcpy failure rather than an error message.
//
// What goes wrong below 4 is NOT the index table: at weight 1, 2 and 3 the table
// holds 5, 17 and 65 entries and every k-mer the seeder can produce indexes
// inside it. A weight-3 run was observed dying on a GPU with
//
//     thrust::system::system_error: device free failed: cudaErrorIllegalAddress
//
// and the mechanism behind that has not been established here. Do not read this
// bound as a fix for that crash; read it as the upstream precondition, enforced.
#define MIN_SEED_WEIGHT 4

// Upper bound on a seed pattern's span (the length of the shape string).
//
// Independent of the weight above: "T" followed by ninety-nine zeroes has a
// legal weight of 1 and a span of 100.  GetKmerIndexAtPos() reads the whole
// span into a fixed nt[MAX_SEED_SPAN] before consulting shape_pos[], so an
// unbounded span writes past that array once per candidate position.
#define MAX_SEED_SPAN 64

#define TRANSITION_MASK 2
#define NUC 8 
#define NUC2 NUC*NUC
#define A_NT 0
#define C_NT 1
#define G_NT 2
#define T_NT 3
#define L_NT 4
#define N_NT 5
#define X_NT 6
#define E_NT 7

#define MAX_BLOCKS 1<<10
#define MAX_THREADS 1024 
#define BLOCK_SIZE 128 
#define NUM_WARPS 4
