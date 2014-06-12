#include <R.h>
#include <Rinternals.h>
#include "genoset.h"

SEXP rangeMeans_rle(SEXP Start, SEXP End, SEXP RunValues, SEXP RunLengths, SEXP Na_rm) {
  int tolerate_na = ! asLogical(Na_rm);
  if (tolerate_na == NA_LOGICAL) { error("'na.rm' must be TRUE or FALSE"); }
  
  int *start_p = INTEGER(Start);
  int *end_p = INTEGER(End);
  int *lengths_p = INTEGER(RunLengths);
  int nrun = (size_t) LENGTH(RunValues);
  int nranges = (size_t) LENGTH(Start);
  
  // Input type dependence
  double *values_p = REAL(RunValues);
  const double na_val = NA_REAL;
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, nranges ));  
  double *ans_p = REAL(Ans);

  // Abstract all the NA checking to a simple lookup of a boolean value
  char* isna = (char *) R_alloc(nrun, sizeof(char));
  isNA(RunValues, isna);
  
  // Just basic C types from here on
  double temp_sum;
  size_t i, start, end, inner_n, sufficient_width, effective_width, run_index;
  size_t lower_run = 0, upper_run = 0;
  size_t* run_start_indices = (size_t*) R_alloc(nrun, sizeof(size_t));
  widthToStart(lengths_p, run_start_indices, nrun);
  size_t last_run = nrun - 1;
  // From here down all type-dependence could be handled by a template on values_p and na_val
  for (i = 0; i < nranges; i++) {
    start = start_p[i];
    end = end_p[i];
    // Find run(s) covered by current range using something like findOverlaps(IRanges(start,width), ranges(rle))
    lower_run = leftBound(run_start_indices, lower_run, last_run, start);
    upper_run = leftBound(run_start_indices, lower_run, last_run, end); // Yes, search the left bound both times
    if (lower_run == upper_run) {  // Range all in one run, special case here allows simpler logic below
      ans_p[i] = values_p[lower_run];
      continue;
    } else {
      // First run
      inner_n = (run_start_indices[lower_run + 1] - start) * !isna[lower_run];
      effective_width = inner_n;
      temp_sum = isna[lower_run] ? 0 : values_p[lower_run] * inner_n;   // floating point NA/Nan contaminate, so have to branch. For na.rm=FALSE, could let them ride. For int types, x * !na would work.
      // Inner runs
      for (run_index = lower_run + 1; run_index < upper_run; run_index++) {
      	inner_n = lengths_p[run_index] * !isna[run_index];
      	effective_width += inner_n;
	temp_sum += isna[run_index] ? 0 : values_p[run_index] * inner_n;   // floating point NA/Nan contaminate, so have to branch. For na.rm=FALSE, could let them ride. For int types, x * !na would work.
      }
      // Last run
      inner_n = ((end - run_start_indices[upper_run]) + 1) * !isna[upper_run];
      effective_width += inner_n;
      temp_sum += isna[upper_run] ? 0 : values_p[upper_run] * inner_n;   // floating point NA/Nan contaminate, so have to branch. For na.rm=FALSE, could let them ride. For int types, x * !na would work.
      // Calculate mean, handling NAs
      sufficient_width = ((start - end) * tolerate_na) + 1;  // Less than this many non-NA values and we return na_val
      ans_p[i] = effective_width < sufficient_width ? na_val : temp_sum / effective_width;
    }
  }
  UNPROTECT(1);
  return Ans;
}

SEXP rangeMeans_numeric(SEXP bounds, SEXP x, SEXP na_rm) {
  if (!isMatrix(bounds) || !isInteger(bounds) || ncols(bounds) != 2) {
    error("'bounds' argument must be a two-column integer matrix.");
  }
  int tolerate_na = ! asLogical(na_rm);
  if (tolerate_na == NA_LOGICAL) { error("'na.rm' must be TRUE or FALSE"); }

  SEXP means, bounds_dimnames, x_dimnames, dimnames;
  int num_cols, num_rows, left, right;

  int num_protected = 0;

  double *x_data = REAL(x);    
  int *bounds_data = INTEGER(bounds);
  int num_bounds = nrows(bounds);
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
  // track left_bound_p, right_bound_p, means_p
  for (int col_index = 0; col_index < num_cols; col_index++) {
    for (int bound_index = 0; bound_index < num_bounds; bound_index++) {
      sum = 0;
      num_na = 0;
      mean_index = (col_index * num_bounds) + bound_index;
      left = bounds_data[bound_index] + col_offset;
      right = bounds_data[bound_index + num_bounds] + col_offset;
      num_to_sum = (right - left) + 1;  // sufficient_width
      for (int i = left-1; i < right; i++) {
	if (! R_FINITE(x_data[i]) ) {  // just NA instead?
	  num_na += 1;
	} else {
	  sum += x_data[i];  // effective_width
	}
      }
      if (num_na == num_to_sum) {
	means_data[ mean_index ] = NA_REAL;
      } else {
	//ans_p[i] = effective_width < sufficient_width ? na_val : temp_sum / effective_width;
	means_data[ mean_index ] = sum / (num_to_sum - num_na);
      }
    }
    col_offset += num_rows;
  }
  
  UNPROTECT(num_protected);
  return(means);
}

SEXP rangeMeans_numeric2(SEXP bounds, SEXP x, SEXP na_rm) {
  if (!isMatrix(bounds) || !isInteger(bounds) || ncols(bounds) != 2) {
    error("'bounds' argument must be a two-column integer matrix.");
  }
  int tolerate_na = ! asLogical(na_rm);
  if (tolerate_na == NA_LOGICAL) { error("'na.rm' must be TRUE or FALSE"); }

  SEXP means, bounds_dimnames, x_dimnames, dimnames;
  int num_cols, num_rows, left, right;
  int num_protected = 0;

  double *x_p = REAL(x);    
  int *bounds_data = INTEGER(bounds);
  int num_bounds = nrows(bounds);
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
  double *means_p = REAL(means);

  double sum;
  int num_na, num_to_sum;
  const double na_val = NA_REAL;
  for (int col_index = 0; col_index < num_cols; col_index++) {
    for (int bound_index = 0; bound_index < num_bounds; means_p++, left_bound++, righ_bound++, bound_index++) {
      sum = 0;
      effective_width = 0;
      left = *left_bound;
      right = *right_bound;
      sufficient_width = ((right - left)*keep_na) + 1;
      for (int i = left-1; i < right; i++) {
	if (! ISNA(x_p[i]) ) {
	  effective_width += 1;
	  sum += x_p[i];
	}
      }
      x_p[i] = effective_width < sufficient_width ? na_val : sum / effective_width;
    }
  }
  
  UNPROTECT(num_protected);
  return(means);
}
