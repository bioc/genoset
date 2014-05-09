#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <R_ext/Utils.h>


// Computing mean of arbitrary slices of an Rle, assuming slices are within [1,length(rle)] (trimmed)
// Takes contents of Rle and IRanges rather than taking these or these stuffed into an RleViews
// The body of this function is nearly independent of special R types, so could perhaps
//   become a separate function, perhaps templated on type of pointer to Values data.
// Having the plain C part a separate function would allow using that function in a loop, 
//  like in a loop over many Rle objects and one IRanges.
SEXP RleViews_viewMeans(SEXP Start, SEXP Width, SEXP Values, SEXP Lengths, SEXP Na_rm) {
  int i, start, width, ans_len, index, lower_run, upper_run, lower_bound, upper_bound, max_index;
  int inner_n, num_na, isna;
  if (!isLogical(Na_rm) || LENGTH(Na_rm) != 1 || LOGICAL(Na_rm)[0] == NA_LOGICAL) {
    error("'na.rm' must be TRUE or FALSE");
  }

  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, LENGTH(Start) ));	
  double *ans_p = REAL(Ans);
  int *start_p = INTEGER(Start);
  int *width_p = INTEGER(Width);
  double *values_p = REAL(Values);
  int *lengths_p = INTEGER(Lengths);
  int keep_na = !LOGICAL(Na_rm)[0];

  // How about abstracting the NA checking to something that just fills out an array of chars with 0,1?  Pass in a char* and the SXP to check, switch on SXP type and loop?
  index = 0;
  upper_run = *lengths_p;
  for (i = 0; i < LENGTH(Start); i++) {
    if (i % 100000 == 99999) { R_CheckUserInterrupt(); }
    start = start_p[i];
    width = width_p[i];
    if (width <= 0) {
      ans_p[i] = R_NaN;
      continue;
    }
    ans_p[i] = 0;
    num_na = 0;
    // I think doing genoset::binary_bound on the cumsums of start and end would be faster here.
    // Computing the inner_n values would then be two subtractions and maybe a loop over any spanned 
    //  runs spanned after the first two.
    //  This would avoid the ternary operator branch (below) most of the time.
    while (index > 0 && upper_run > start) {
      upper_run -= *lengths_p;
      lengths_p--;
      index--;
    }
    while (upper_run < start) {
      lengths_p++;
      index++;
      upper_run += *lengths_p;
    }
    lower_run = upper_run - *lengths_p + 1;
    lower_bound = start;
    upper_bound = start + width - 1;
    while (lower_run <= upper_bound) {
      isna = ISNA(values_p[index]); // Depends on TYPEOF(Values), abstract to char array of 0,1
      inner_n = (1 + (upper_bound < upper_run ? upper_bound : upper_run) - (lower_bound > lower_run ? lower_bound : lower_run));
      num_na += isna * inner_n;
      ans_p[i] += values_p[index] * (inner_n * (!isna)); // Depends on typeof(values_p), polymorphism with template?
      if (index < max_index) {
	lengths_p++;
	index++;
	lower_run = upper_run + 1;
	lower_bound = lower_run;
	upper_run += *lengths_p;
      } else {
	break;
      }
    }
    if ( num_na == width || (keep_na && num_na > 0)) {
      ans_p[i] = NA_REAL;
    } else {
      ans_p[i] /= (width - num_na);
    }
  }
  UNPROTECT(1);
  return Ans;
}
