#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <R_ext/Utils.h>

SEXP RleViews_viewMeans(SEXP Start, SEXP Width, SEXP Values, SEXP Lengths, SEXP Na_rm) {
  int i, start, width, ans_len, index, lower_run, upper_run, lower_bound, upper_bound, max_index;
  int inner_n, num_na, isna;
  if (!isLogical(Na_rm) || LENGTH(Na_rm) != 1 || LOGICAL(Na_rm)[0] == NA_LOGICAL) {
    error("'na.rm' must be TRUE or FALSE");
  }
  
  int keep_na = !LOGICAL(Na_rm)[0];
  ans_len = LENGTH(Start);
  SEXP Ans;
  PROTECT(Ans = allocVector(REALSXP, ans_len));	
  double *ans_p = REAL(Ans);
  int *start_p = INTEGER(Start);
  int *width_p = INTEGER(Width);
  int *lengths_p = INTEGER(Lengths);
  // How about abstracting the NA checking to something that just fills out an array of chars with 0,1?  Pass in a char* and the SXP to check, switch on SXP type and loop?
  index = 0;
  upper_run = *lengths_p;
  for (i = 0; i < ans_len; i++) {
    if (i % 100000 == 99999) { R_CheckUserInterrupt(); }
    start = start_p[i];
    width = width_p[i];
    if (width <= 0) {
      ans_p[i] = R_NaN;
      continue;
    }
    ans_p[i] = 0;
    num_na = 0;
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
      isna = ISNA(REAL(Values)[index]);
      inner_n = (1 + (upper_bound < upper_run ? upper_bound : upper_run) - (lower_bound > lower_run ? lower_bound : lower_run));
      num_na += isna * inner_n;
      ans_p[i] += REAL(Values)[index] * (inner_n * (!isna));
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



