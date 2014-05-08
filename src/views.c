#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <R_ext/Utils.h>

SEXP RleViews_viewMeans(SEXP Start, SEXP Width, SEXP Values, SEXP Lengths, SEXP Na_rm) {

	int i, n, start, width, ans_len, index, lower_run, upper_run, lower_bound, upper_bound, max_index;

	if (!isLogical(Na_rm) || LENGTH(Na_rm) != 1 || LOGICAL(Na_rm)[0] == NA_LOGICAL) {
	  error("'na.rm' must be TRUE or FALSE");
	}
	
	int do_na = !LOGICAL(Na_rm)[0];
	ans_len = LENGTH(Start);
	SEXP Ans;
	PROTECT(Ans = allocVector(REALSXP, ans_len));	
	double *ans_p = REAL(Ans);
	int *start_p = INTEGER(Start);
	int *width_p = INTEGER(Width);
	int *lengths_p = INTEGER(Lengths);
//	int ans_is_real = 0;
//	switch (TYPEOF(Values)) {
//	  case LGLSXP:
//	  case INTSXP:
//	    ans_is_real = 0;
//	    break;
//	  case REALSXP:
//	    ans_is_real = 1;
//	    break;
//	  default:
//	    error("Rle must contain either 'integer' or 'numeric' values");
//	  }
	index = 0;
	upper_run = *lengths_p;
	for (i = 0; i < ans_len; i++) {
	  if (i % 100000 == 99999) { R_CheckUserInterrupt(); }
	  start = start_p[i];
	  width = width_p[i];
	  //	  printf("Width: %i\n",width);
	  if (width <= 0) {
	    ans_p[i] = R_NaN;
	  } else {
	    n = width;
	    ans_p[i] = 0;
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
	      //	      if (ans_is_real ? ISNA(REAL(Values)[index]) : INTEGER(Values)[index] == NA_INTEGER) {
	      if (ISNA(REAL(Values)[index])) {
		if (do_na) {
		  ans_p[i] = NA_REAL;
		  break;
		}
		n -= (1 + (upper_bound < upper_run ? upper_bound : upper_run) -
		   (lower_bound > lower_run ? lower_bound : lower_run));
		
	      } else {
		//		  ans_p[i] += (ans_is_real ? REAL(Values)[index] : INTEGER(Values)[index]) *
		  ans_p[i] += REAL(Values)[index] *
		  (1 + (upper_bound < upper_run ? upper_bound : upper_run) -
		   (lower_bound > lower_run ? lower_bound : lower_run));
	      }
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
	    if (n == 0) {
	      ans_p[i] = R_NaN;
	    } else if (ans_p[i] != NA_REAL) {
	      ans_p[i] /= n;
	    }
	  }
	}
	UNPROTECT(1);
	return Ans;
}
