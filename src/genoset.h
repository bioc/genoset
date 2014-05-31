#ifndef GENOSET_H
#define GENOSET_H

// utils.c
void isNA(SEXP vec, char* na);
int numNA(SEXP vec, char* na); // Would rather overload (C++ only) with int isNA(SEXP vec, char* na);
void widthToStartEnd(int* width, double* start, double* end, int n);
void widthToStart(int* width, int* start, int n);
int leftBound(int* positions, int query, int n, unsigned int restart);

// views.c
SEXP RleViews_viewMeans(SEXP Start, SEXP Width, SEXP Values, SEXP Lengths, SEXP Na_rm);

// bounds.c
SEXP binary_bound( SEXP starts, SEXP stops, SEXP positions, SEXP valid_indices);
SEXP rangeColMeans( SEXP bounds, SEXP x );
SEXP binary_bound_by_chr(SEXP nquery, SEXP query_chr_indices, SEXP query_starts, SEXP query_ends, SEXP query_names, SEXP subject_chr_indices, SEXP subject_starts, SEXP subject_ends);

#endif

