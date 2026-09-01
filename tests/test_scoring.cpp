// Tests that a LASTZ scoring file's bad_score and fill_score reach the
// substitution matrix the GPU scores against.
//
// Run with:  bash tests/test_scoring.bash
//
// Needs a C++ compiler, but no CUDA toolkit and no GPU: it drives the real
// LASTZ score-set reader in common/dna_utilities.c through common/scoring.c.
//
// ---------------------------------------------------------------------------
// The bug this exists for: --scoring accepted a LASTZ scoring file but used
// only its 4x4 ACGT block. bad_score and fill_score were parsed by the LASTZ
// reader and then discarded, and hardcoded -1000/-100 used instead -- while
// runner.py passes the SAME file to lastz via --scores=, where LASTZ does
// honour them. One file, two interpretations, inside one pipeline.
// ---------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>

#include "scoring.h"

static int failures = 0;

static void expect_eq(const char* label, int actual, int expected)
{
    if (actual != expected) {
        printf("not ok - %s\n     got: %d\nexpected: %d\n", label, actual, expected);
        failures++;
    } else {
        printf("ok - %s\n", label);
    }
}

int main(int argc, char** argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s <scoring-file>\n", argv[0]); return 2; }

    const int xdrop = 910;

    // --- the default path, with no scoring file -----------------------------
    // Pins that moving this out of main() changed nothing.
    int hoxd70[L_NT][L_NT] = {{   91, -114,  -31, -123},
                              { -114,  100, -125,  -31},
                              {  -31, -125,  100, -114},
                              { -123,  -31, -114,   91}};
    int sub_mat[NUC*NUC];
    build_substitution_matrix(sub_mat, hoxd70, -1000, -100, "x", -100, -100, xdrop);
    expect_eq("default: A/A keeps the HOXD70 score",   sub_mat[A_NT*NUC+A_NT],  91);
    expect_eq("default: lower case is bad_score",      sub_mat[A_NT*NUC+L_NT], -1000);
    expect_eq("default: N is bad_score",               sub_mat[N_NT*NUC+A_NT], -1000);
    expect_eq("default: other IUPAC is fill_score",    sub_mat[A_NT*NUC+X_NT], -100);
    expect_eq("default: separator is -10*xdrop",       sub_mat[A_NT*NUC+E_NT], -9100);

    // --- --ambiguous=n ------------------------------------------------------
    build_substitution_matrix(sub_mat, hoxd70, -1000, -100, "n", 0, 0, xdrop);
    expect_eq("ambiguous=n: N substitutions score 0",  sub_mat[N_NT*NUC+A_NT], 0);
    expect_eq("ambiguous=n: N/N takes the reward",     sub_mat[N_NT*NUC+N_NT], 0);

    // --- a scoring file carrying its own bad_score and fill_score ------------
    int from_file[L_NT][L_NT];
    int bad_score = 0, fill_score = 0;
    load_scoring_matrix(from_file, &bad_score, &fill_score, argv[1]);

    expect_eq("scoring file: A/A read from the file",  from_file[A_NT][A_NT],  85);
    expect_eq("scoring file: A/C read from the file",  from_file[A_NT][C_NT], -110);
    expect_eq("scoring file: bad_score read from the file",  bad_score,  -750);
    expect_eq("scoring file: fill_score read from the file", fill_score,  -55);

    build_substitution_matrix(sub_mat, from_file, bad_score, fill_score, "x", -100, -100, xdrop);
    expect_eq("scoring file: lower case uses the file's bad_score",  sub_mat[A_NT*NUC+L_NT], -750);
    expect_eq("scoring file: N uses the file's bad_score",           sub_mat[N_NT*NUC+A_NT], -750);
    expect_eq("scoring file: other IUPAC uses the file's fill_score",sub_mat[A_NT*NUC+X_NT],  -55);
    expect_eq("scoring file: X/X uses the file's fill_score",        sub_mat[X_NT*NUC+X_NT],  -55);

    if (failures > 0) { printf("\n%d test(s) failed\n", failures); return 1; }
    printf("\nall tests passed\n");
    return 0;
}
