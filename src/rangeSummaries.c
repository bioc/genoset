 #include <R.h>
 #include <Rinternals.h>
 #include "genoset.h"

 // Computing mean of arbitrary slices of an Rle, assuming slices are within [1,length(rle)] (trimmed)
 // Takes contents of Rle and IRanges rather than taking these or these stuffed into an RleViews
 // The body of this function is nearly independent of special R types, so could perhaps
 //   become a separate function, perhaps templated on type of pointer to Values data.
 // Having the plain C part a separate function would allow using that function in a loop, 
 //  like in a loop over many Rle objects and one IRanges.
 // Assuming Values gets as.numeric on the way in. No support for complex type.  For mean, we have to cast to double one at a time anyway.

// Some branching for NA check only necessary for na.rm=TRUE and numeric Values.
// Could let the NaNs contaminate for na.rm=FALSE or use x * !na[i] for int types, 
//   but this is a minor fraction of the time.
SEXP rangeMeans_rle(SEXP Start, SEXP Width, SEXP Values, SEXP Lengths, SEXP Na_rm) {
  int keep_na = ! asLogical(Na_rm);
  if (keep_na == NA_LOGICAL) { error("'na.rm' must be TRUE or FALSE"); }
  
  int *start_p = INTEGER(Start);
  int *width_p = INTEGER(Width);
  int *lengths_p = INTEGER(Lengths);
  int nrun = LENGTH(Values);
  int nranges = LENGTH(Start);
  
  // Input type dependence
  double *values_p = REAL(Values);
  const double na_val = NA_REAL;
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, nranges ));  
  double *ans_p = REAL(Ans);

  // Abstract all the NA checking to a simple lookup of a boolean value
  char* isna = (char *) R_alloc(nrun, sizeof(char));
  isNA(Values, isna);
  
  // Just basic C types from here on
  double temp_sum;
  int i, start, width, end, inner_n, effective_width, run_index;
  int lower_run = 0, upper_run = 0;
  int* run_first_index = (int*) R_alloc(nrun, sizeof(int));
  widthToStart(lengths_p, run_first_index, nrun);
  // From here down all type-dependence could be handled by a template on values_p and na_val
  for (i = 0; i < nranges; i++) {
    start = start_p[i];
    width = width_p[i];
    end = (start + width) - 1;
    // Find run(s) covered by current range using something like findOverlaps(IRanges(start,width), ranges(rle))
    lower_run = leftBound(run_first_index, start, nrun, lower_run);
    upper_run = leftBound(run_first_index, end, nrun, lower_run);  // Yes, search the left bound both times
    //    printf("lower_run: %i, upper_run: %i\n", lower_run, upper_run);
    if (lower_run == upper_run) {  // Range all in one run, special case here allows simpler logic below
      ans_p[i] = values_p[lower_run];
      continue;
    } else {
      // First run
      inner_n = (run_first_index[lower_run + 1] - start) * !isna[lower_run];
      effective_width = inner_n;
      temp_sum = isna[lower_run] ? 0 : values_p[lower_run] * inner_n;   // floating point NA/Nan contaminate, so have to branch. For na.rm=FALSE, could let them ride. For int types, x * !na would work.
      // Inner runs
      for (run_index = lower_run + 1; run_index < upper_run; run_index++) {
      	inner_n = lengths_p[run_index] * !isna[run_index];
      	effective_width += inner_n;
	temp_sum += isna[run_index] ? 0 : values_p[run_index] * inner_n;   // floating point NA/Nan contaminate, so have to branch. For na.rm=FALSE, could let them ride. For int types, x * !na would work.
      }
      // Last run
      inner_n = ((end - run_first_index[upper_run]) + 1) * !isna[upper_run];
      effective_width += inner_n;
      temp_sum += isna[upper_run] ? 0 : values_p[upper_run] * inner_n;   // floating point NA/Nan contaminate, so have to branch. For na.rm=FALSE, could let them ride. For int types, x * !na would work.
      ans_p[i] = (effective_width != width && (effective_width == 0 || keep_na)) ? na_val : temp_sum / effective_width;
    }
  }

  UNPROTECT(1);
  return Ans;
}

SEXP rangeMeans_vector( SEXP bounds, SEXP x ) {
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
  int num_na, num_to_sum;
  int col_offset = 0;
  int mean_index = 0;
  for (int col_index = 0; col_index < num_cols; col_index++) {
    for (int bound_index = 0; bound_index < num_bounds; bound_index++) {
      sum = 0;
      num_na = 0;
      mean_index = (col_index * num_bounds) + bound_index;
      left = bounds_data[bound_index] + col_offset;
      right = bounds_data[bound_index + num_bounds] + col_offset;
      num_to_sum = (right - left) + 1;
      for (int i = left-1; i < right; i++) {
  if (! R_FINITE(x_data[i]) ) {
	  num_na += 1;
	} else {
	  sum += x_data[i];
	}
      }
      if (num_na == num_to_sum) {
	means_data[ mean_index ] = NA_REAL;
      } else {
	means_data[ mean_index ] = sum / (num_to_sum - num_na);
      }
    }
    col_offset += num_rows;
  }
  
  UNPROTECT(num_protected);
  return(means);
}
