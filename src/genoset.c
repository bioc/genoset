#include <R.h>
#include <Rinternals.h>

void binary_bound( int *starts, int *stops, int *positions, int *num_queries, int *num_positions, int *bounds, int *valid_indices, int *offset ) {
  int query_index, probe, left, right, low, high, jump;
  low = -1; /*  Set low off left end */
  high = num_positions[0]; /* Set high to off right end */
  
  for (query_index=0; query_index < num_queries[0]; query_index++) {
   
    /* Left bound */
    left = starts[query_index];

    /* If data unsorted, current target may be anywhere left of low for previous target, just start at left */
    if (low >= 0 && left < positions[low]) {
      low = -1;
    }
    
    /* Right bound likely close to high from previous gene */
    for (jump=1; ; jump*=2) {
      high += jump;
      if (high >= num_positions[0]) {
	high = num_positions[0];
	break;
      }
      if (left < positions[high]) {  /* Note difference to similar code resetting high after first binary search */
	break;
      }
      low = high;
    }

    /* Now binary search for right bound */
    while (high - low > 1) {
      probe = (high + low) >> 1;
      if (positions[probe] > left) {
        high = probe;
      } else {
	low = probe;
      }
    }
    bounds[query_index] = low;
    
    /* Right bound */
    right = stops[query_index];

    /* Right bound likely close to left bound, relative to length of positions, so expand exponentially, a la findInterval */
    for (jump=1; ; jump*=2) {
      high += jump;
      if (high >= num_positions[0]) {
	high = num_positions[0];
	break;
      }
      if (right <= positions[high]) {
	break;
      }
      low = high;
    }

    /* Now binary search for right bound */
    while (high - low > 1) {
      probe = (high + low) >> 1;
      if (positions[probe] < right) {
        low = probe;
      } else {
        high = probe;
      }
    }
    bounds[query_index + num_queries[0]] = high; /* Right bound goes in second column of bounds */
  
    low = bounds[query_index]; /* Reset low to left end of this query to start next query */
    
  } /* End foreach query */
  
  offset[0] += 1;
  for(int i=0; i < num_queries[0]; i++) {
    if (valid_indices[0] == 1) {
      if (bounds[i] < 0) {
	bounds[i] = 0;
      }
      if (bounds[i + num_queries[0]] >= num_positions[0]) {
	bounds[i + num_queries[0]] = num_positions[0] - 1;
      }
    }
    bounds[i] += offset[0];
    bounds[i + num_queries[0]] += offset[0];
  }
}

SEXP rangeColMeans( SEXP bounds, SEXP x ) {
  SEXP means, bounds_dimnames, x_dimnames, dimnames;
  int num_cols, num_rows, left, right;

  int num_protected = 0;

  double *x_data = REAL(x);    
  int *bounds_data = INTEGER(bounds);
  int num_bounds = length(bounds) / 2;
  bounds_dimnames = getAttrib(bounds, R_DimNamesSymbol);
  x_dimnames = getAttrib(x, R_DimNamesSymbol);

  if (isMatrix(x)) {
    num_cols = ncols(x);
    num_rows = nrows(x);
    PROTECT(means = allocMatrix(REALSXP, num_bounds, num_cols)); num_protected++;
    if ( GetRowNames(bounds_dimnames) != R_NilValue || GetColNames(x_dimnames) != R_NilValue) {
    PROTECT(dimnames = allocVector(VECSXP, 2)); num_protected++;
    SET_VECTOR_ELT(dimnames, 0, duplicate(GetRowNames(bounds_dimnames)));
    SET_VECTOR_ELT(dimnames, 1, duplicate(GetColNames(x_dimnames)));
    setAttrib(means, R_DimNamesSymbol, dimnames);
    }
  } else {
    num_cols = 1;
    num_rows = length(x);
    PROTECT(means = allocVector(REALSXP, num_bounds)); num_protected++;
    if ( GetRowNames(bounds_dimnames) != R_NilValue ) {
      setAttrib(means, R_NamesSymbol, duplicate(GetRowNames(bounds_dimnames)));
    }
  }
  double *means_data = REAL(means);

  double sum;
  int col_offset = 0;
  int mean_index = 0;
  for (int col_index = 0; col_index < num_cols; col_index++) {
    for (int bound_index = 0; bound_index < num_bounds; bound_index++) {
      sum = 0;
      mean_index = (col_index * num_bounds) + bound_index;
      left = bounds_data[bound_index] + col_offset;
      right = bounds_data[bound_index + num_bounds] + col_offset;
      for (int i = left-1; i < right; i++) {
	sum += x_data[i];
      }
      means_data[ mean_index ] = sum / ((right - left)+1);
    }
    col_offset += num_rows;
  }
  
  UNPROTECT(num_protected);
  return(means);
}
