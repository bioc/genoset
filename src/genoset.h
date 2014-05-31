#ifndef GENOSET_H
#define GENOSET_H

// utils.c
void isNA(SEXP vec, char* na);
int numNA(SEXP vec, char* na); // Would rather overload (C++ only) with int isNA(SEXP vec, char* na);
void widthToStartEnd(int* width, int* start, int* end, int n);
void widthToStart(int* width, int* start, int n);
int leftBound(int* positions, int query, int n, unsigned int restart);

// views.c
SEXP rangeMeans_rle(SEXP Start, SEXP Width, SEXP Values, SEXP Lengths, SEXP Na_rm);
SEXP rangeMeans_vector(SEXP bounds, SEXP x);

// bounds.c
SEXP binary_bound( SEXP starts, SEXP stops, SEXP positions, SEXP valid_indices);
SEXP binary_bound_by_chr(SEXP nquery, SEXP query_chr_indices, SEXP query_starts, SEXP query_ends, SEXP query_names, SEXP subject_chr_indices, SEXP subject_starts, SEXP subject_ends);

#endif

