#include <string.h>

#include "scoring.h"

// LASTZ's own defaults, from the score set file format documented in
// dna_utilities.c: bad_score is the score for the row and column of the "bad"
// character (X, for DNA), fill_score is used for any pair the matrix does not
// define.
#define DEFAULT_BAD_SCORE  -1000
#define DEFAULT_FILL_SCORE  -100

static int in_uchars(const u8* s, int ch) {
	for (; *s != 0; s++) {
		if ((int) *s == ch) {
			return 1;
		}
	}
	return 0;
}

// LASTZ stores bad_score by writing it across the row and column of the bad
// character (X, for DNA), rather than keeping the scalar.
static int recover_bad_score(scoreset* ss) {
	if (ss->badRow == 0 || ss->badCol == 0) {
		return DEFAULT_BAD_SCORE;
	}
	return (int) ss->sub[ss->badRow][ss->badCol];
}

// LASTZ fills all 256x256 entries with fill_score and then writes the matrix,
// the lower case copies and the bad row/column over the top, so any pair it
// never touched still holds fill_score.
static int recover_fill_score(scoreset* ss) {
	int r, c;

	for (r = 1; r < 256; r++) {
		if (r == ss->badRow || in_uchars(ss->rowChars, r)) {
			continue;
		}
		for (c = 1; c < 256; c++) {
			if (c == ss->badCol || in_uchars(ss->colChars, c)) {
				continue;
			}
			return (int) ss->sub[r][c];
		}
	}

	return DEFAULT_FILL_SCORE;
}

void load_scoring_matrix(int scoring_matrix[][L_NT], int* bad_score, int* fill_score, char* scoreFilename) {
	exscoreset*     xss = NULL;
	u8*             rowChars = (u8*)"ACGT";
	u8*             colChars = (u8*)"ACGT";
	u8*             r, *c;
	int             x, y;

	*bad_score  = DEFAULT_BAD_SCORE;
	*fill_score = DEFAULT_FILL_SCORE;

	xss = read_score_set_by_name (scoreFilename);

	for (r = rowChars, x = 0; *r != 0; r++, x++) {
		for (c = colChars, y = 0; *c != 0; c++, y++) {
			scoring_matrix[x][y]  = (int) ((scoreset*) xss)->sub[*r][*c];
		}
	}

	*bad_score  = recover_bad_score((scoreset*) xss);
	*fill_score = recover_fill_score((scoreset*) xss);

	// No NULL check: read_score_set_by_name() calls suicide()/fopen_or_die() on
	// every failure path, so it either returns a score set or does not return.
	free_if_valid("score set", xss);
}

void build_substitution_matrix(int* sub_mat, int scoring_matrix[][L_NT], int bad_score, int fill_score, const char* ambiguous_field, int ambiguous_reward, int ambiguous_penalty, int xdrop) {
	int i, j;

	//ACGT
	for(i = 0; i < L_NT; i++){
		for(j = 0; j < L_NT; j++){
			sub_mat[i*NUC+j] = scoring_matrix[i][j];
		}
	}

	//lower case characters
	for(i = 0; i < L_NT; i++){
		sub_mat[i*NUC+L_NT] = bad_score;
		sub_mat[L_NT*NUC+i] = bad_score;
	}
	sub_mat[L_NT*NUC+L_NT] = bad_score;

	//N
	if(strcmp(ambiguous_field, "n") == 0 || strcmp(ambiguous_field, "iupac") == 0){
		for(i = 0; i < N_NT; i++){
			sub_mat[i*NUC+N_NT] = ambiguous_penalty;
			sub_mat[N_NT*NUC+i] = ambiguous_penalty;
		}
		sub_mat[N_NT*NUC+N_NT] = ambiguous_reward;
	}
	else{
		for(i = 0; i < N_NT; i++){
			sub_mat[i*NUC+N_NT] = bad_score;
			sub_mat[N_NT*NUC+i] = bad_score;
		}
		sub_mat[N_NT*NUC+N_NT] = bad_score;
	}

	//other IUPAC
	if(strcmp(ambiguous_field, "iupac") == 0){
		for(i = 0; i < X_NT; i++){
			sub_mat[i*NUC+X_NT] = ambiguous_penalty;
			sub_mat[X_NT*NUC+i] = ambiguous_penalty;
		}
		sub_mat[X_NT*NUC+X_NT] = ambiguous_reward;
	}
	else{
		for(i = 0; i < L_NT; i++){
			sub_mat[i*NUC+X_NT] = fill_score;
			sub_mat[X_NT*NUC+i] = fill_score;
		}

		for(i = L_NT; i < X_NT; i++){
			sub_mat[i*NUC+X_NT] = bad_score;
			sub_mat[X_NT*NUC+i] = bad_score;
		}
		sub_mat[X_NT*NUC+X_NT] = fill_score;
	}

	//block separator
	for(i = 0; i < E_NT; i++){
		sub_mat[i*NUC+E_NT] = -10*xdrop;
		sub_mat[E_NT*NUC+i] = -10*xdrop;
	}
	sub_mat[E_NT*NUC+E_NT] = -10*xdrop;
}
