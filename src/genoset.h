#ifndef GENOSET_H
#define GENOSET_H

// utils.c
void isNA(SEXP vec, char* na);
int numNA(SEXP vec, char* na); // Would rather overload (C++ only) with int isNA(SEXP vec, char* na);
void widthToStartEnd(int* width, double* start, double* end, int n);
void widthToStart(int* width, double* start, int n);

// views.c
SEXP RleViews_viewMeans(SEXP Start, SEXP Width, SEXP Values, SEXP Lengths, SEXP Na_rm);

// genoset.c
void binary_bound( int *starts, int *stops, int *positions, int *num_queries, int *num_positions, int *bounds, int *valid_indices, int *offset );
SEXP rangeColMeans( SEXP bounds, SEXP x );
SEXP binary_bound_by_chr(SEXP nquery, SEXP query_chr_indices, SEXP query_starts, SEXP query_ends, SEXP query_names, SEXP subject_chr_indices, SEXP subject_starts, SEXP subject_ends);

#endif

