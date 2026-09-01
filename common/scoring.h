#ifndef scoring_H                             // (prevent multiple inclusion)
#define scoring_H

#ifdef __cplusplus
extern "C" {
#endif

#include "dna_utilities.h"
#include "parameters.h"

// Read the ACGT substitution matrix from a LASTZ scoring file, along with the
// bad_score and fill_score settings the same file may carry.
void load_scoring_matrix(int scoring_matrix[][L_NT], int* bad_score, int* fill_score, char* scoreFilename);

// Expand the 4x4 ACGT matrix into the NUC x NUC table the GPU scores against,
// filling in the rows for lower case, N, other IUPAC codes and the block
// separator.  ambiguous_field is "x", "n" or "iupac" (the first field of
// --ambiguous).
void build_substitution_matrix(int* sub_mat, int scoring_matrix[][L_NT], int bad_score, int fill_score, const char* ambiguous_field, int ambiguous_reward, int ambiguous_penalty, int xdrop);

#ifdef __cplusplus
}
#endif

#endif // scoring_H
