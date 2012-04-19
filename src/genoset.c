#include <R.h>
#include <Rinternals.h>

/* TODO: break out redundant code in binary_bound and binary_bound_by_chr for ease of reading/maintenance */

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
  int num_na;
  int col_offset = 0;
  int mean_index = 0;
  for (int col_index = 0; col_index < num_cols; col_index++) {
    for (int bound_index = 0; bound_index < num_bounds; bound_index++) {
      sum = 0;
      num_na = 0;
      mean_index = (col_index * num_bounds) + bound_index;
      left = bounds_data[bound_index] + col_offset;
      right = bounds_data[bound_index + num_bounds] + col_offset;
      for (int i = left-1; i < right; i++) {
	if (! R_FINITE(x_data[i]) ) {
	  num_na += 1;
	} else {
	  sum += x_data[i];
	}
      }
      means_data[ mean_index ] = sum / (((right - left)+1)-num_na);
    }
    col_offset += num_rows;
  }
  
  UNPROTECT(num_protected);
  return(means);
}

SEXP binary_bound_by_chr(SEXP nquery, SEXP query_chr_indices, SEXP query_starts, SEXP query_ends, SEXP query_names, SEXP subject_chr_indices, SEXP subject_starts, SEXP subject_ends) {

  /* Pointers into data in arg objects */
  double *p_query_chr_indices = REAL(query_chr_indices);
  double *p_subject_chr_indices = REAL(subject_chr_indices);
  int *p_query_starts = INTEGER(query_starts);
  int *p_query_ends = INTEGER(query_ends);
  int *p_subject_starts = INTEGER(subject_starts);
  int *p_subject_ends = INTEGER(subject_ends);
  p_query_starts--; p_query_ends--;  // -- gives 1-based indexing
  p_subject_starts--; p_subject_ends--;  // -- gives 1-based indexing

  /* Make results matrix, "bounds" */
  int num_protected = 0;
  SEXP bounds, dimnames, rownames, colnames;
  PROTECT(bounds = allocMatrix(INTSXP, INTEGER(nquery)[0], 2)); num_protected++;
  PROTECT(dimnames = allocVector(VECSXP, 2)); num_protected++;
  PROTECT(rownames = allocVector(STRSXP, INTEGER(nquery)[0])); num_protected++;
  PROTECT(colnames = allocVector(STRSXP, 2)); num_protected++;
  SET_VECTOR_ELT(dimnames, 0, rownames);
  SET_STRING_ELT(colnames, 0, mkChar("left"));
  SET_STRING_ELT(colnames, 1, mkChar("right"));
  SET_VECTOR_ELT(dimnames, 1, colnames);
  setAttrib(bounds, R_DimNamesSymbol, dimnames);
  int *p_bounds = INTEGER(bounds);

  int low, high, probe, jump;  /* markers for binary search */
  int init_low, init_high;  /* index bounds for search in chr */
  int left, right;  /* values to search for */

  /* markers for row in two cols of result matrix, bounds */
  int left_bound_index = 0;
  int right_bound_index = INTEGER(nquery)[0];

  /* foreach chr */
  for (int chr_start_index=0, chr_end_index=nrows(query_chr_indices); chr_start_index < nrows(query_chr_indices); chr_start_index++, chr_end_index++) {
    /* binary_bound all queries in chr*/
    /* How best to handle being off each end of chr?  Test upfront?  Test at end? */

    init_low = low = p_subject_chr_indices[chr_start_index];
    init_high = high = p_subject_chr_indices[chr_end_index];

    /* This loop should work on query_index in range specified by query_chr_indices for this chr, left and right bound indices must increment every loop but be set at the top */
    for (int query_index=p_query_chr_indices[chr_start_index]; query_index <= p_query_chr_indices[chr_end_index]; query_index++, left_bound_index++, right_bound_index++) {
      SET_STRING_ELT(rownames,left_bound_index, STRING_ELT(query_names,query_index-1));
      left = p_query_starts[query_index];
      right = p_query_ends[query_index];

      // Check for being completely off either end
      if (right <= p_subject_starts[init_low]) {
	p_bounds[right_bound_index] = init_low;
	p_bounds[left_bound_index] = init_low;
	continue;
      }
      if (left >= p_subject_ends[init_high]) {
	p_bounds[right_bound_index] = init_high;
	p_bounds[left_bound_index] = init_high;
	continue;
      }

      /* Left bound */
      /* If data unsorted, current target may be anywhere left of low for previous target, just start at left */
      if (left < p_subject_starts[low]) {
	low = init_low;
      }

      /* Right bound likely close to high from previous gene */
      for (jump=1; ; jump*=2) {
	high += jump;
	if (high >= init_high) {
	  high = init_high;
	  break;
	}
	if (left < p_subject_starts[high]) {  /* Note difference to similar code resetting high after first binary search */
	  break;
	}
	low = high;
      }

      /* Now binary search for left bound */
      while (high - low > 1) {
	probe = (high + low) >> 1;

	if (p_subject_starts[probe] > left) {
	  high = probe;
	} else {
	  low = probe;
	}
      }
      p_bounds[left_bound_index] = low;

      /* Right bound */
      /* Right bound likely close to left bound, relative to length of positions, so expand exponentially, a la findInterval */
      for (jump=1; ; jump*=2) {
	high += jump;
	if (high >= init_high) {
	  high = init_high;
	  break;
	}
	if (right <= p_subject_ends[high]) {
	  break;
	}
	low = high;
      }

      /* Now binary search for right bound */
      while (high - low > 1) {
	probe = (high + low) >> 1;
	if (p_subject_ends[probe] < right) {
	  low = probe;
	} else {
	  high = probe;
	}
      }
      p_bounds[right_bound_index] = high; /* Right bound goes in second column of bounds */
  
      low = p_bounds[left_bound_index]; /* Reset low to left end of this query to start next query */
    
    } /* End foreach query in chr*/
  } /* End foreach chr */
  UNPROTECT(num_protected);
  return(bounds);
}
