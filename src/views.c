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
 SEXP RleViews_viewMeans(SEXP Start, SEXP Width, SEXP Values, SEXP Lengths, SEXP Na_rm) {
   if (!isLogical(Na_rm) || LENGTH(Na_rm) != 1 || LOGICAL(Na_rm)[0] == NA_LOGICAL) {
     error("'na.rm' must be TRUE or FALSE");
   }

   int length_start = LENGTH(Start);
   SEXP Ans;
   PROTECT(Ans = allocVector(REALSXP, length_start ));	
   double *ans_p = REAL(Ans);
   int *start_p = INTEGER(Start);
   int *width_p = INTEGER(Width);
   double *values_p = REAL(Values);
   int *lengths_p = INTEGER(Lengths);
   int keep_na = !LOGICAL(Na_rm)[0];

   double temp_sum;
   int i, start, width, ans_len, index, lower_run, upper_run, lower_bound, upper_bound, max_index;
   int inner_n, num_na, isna;
   // How about abstracting the NA checking to something that just fills out an array of chars with 0,1?  Pass in a char* and the SXP to check, switch on SXP type and loop?
   max_index = LENGTH(Lengths) - 1;
   index = 0;
   upper_run = *lengths_p;
   for (i = 0; i < length_start; i++) {
     start = start_p[i];
     width = width_p[i];
     if (width <= 0) {
       ans_p[i] = R_NaN;
       continue;
     }
     temp_sum = 0;
     num_na = 0;
     // I think doing genoset::binary_bound on the cumsums of start and end would be faster here for finding lowe and upper run.
     // Hmm, looks like bound finding is linear for first, but then looks nearby for subsequent.  Maybe not much gain with binary search with sorted ranges.
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
       isna = ISNA(values_p[index]); // Depends on TYPEOF(Values), abstract to char array of 0,1 with genoset/src/utils.c isNA
       inner_n = 1 + (upper_bound - lower_bound);
       num_na += isna * inner_n;
       temp_sum += values_p[index] * (inner_n * (!isna)); // Depends on typeof(temp_sum) or typeof(values_p)
       if (index >= max_index) { break; }
       lengths_p++;
       index++;
       lower_run = upper_run + 1;
       lower_bound = lower_run;
       upper_run += *lengths_p;
     }
     if ( num_na == width || (keep_na && num_na > 0)) {
       temp_sum = NA_REAL;
     } else {
       temp_sum /= (width - num_na);
     }
     ans_p[i] = temp_sum;
   }
   UNPROTECT(1);
   return Ans;
 }

SEXP RleViews_viewMeans2(SEXP Start, SEXP Width, SEXP Values, SEXP Lengths, SEXP Na_rm) {
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
  
  double temp_sum;
  int i, start, width, inner_n, effective_width, isna;
  int lower_run, upper_run, run_index, mflag = 0;
  
  double* run_first_index = (double *) R_alloc(nrun, sizeof(double));
  double* run_last_index = (double *) R_alloc(nrun, sizeof(double));

  widthToStartEnd(lengths_p, run_first_index, run_last_index, nrun);
  
  for (i = 0; i < nranges; i++) {
    start = start_p[i];
    width = width_p[i];
    temp_sum = 0;
    effective_width = 0;
    // Find run(s) covered by current range using something like findOverlaps(IRanges(start,width), ranges(rle))
    lower_run = findInterval(run_first_index, nrun, start, 0, 0, lower_run, &mflag) - 1;
    upper_run = findInterval(run_first_index, nrun, (start + width) - 1, 0, 0, lower_run, &mflag) - 1;  // Yes, search the left bound both times
    if (lower_run == upper_run) {  // Range all in one run, special case here allows simpler logic below
      ans_p[i] = values_p[lower_run];
      continue;
    } else {
      // First run
      isna = ISNA(values_p[lower_run]); // Depends on TYPEOF(Values), abstract to char array of 0,1 
      inner_n = (1 + (run_last_index[lower_run] - start)) * !isna;
      effective_width += inner_n;
      temp_sum += values_p[lower_run] * inner_n;
      // Inner runs
      for (run_index = lower_run + 1; run_index < upper_run; run_index++) {
      	isna = ISNA(values_p[run_index]);
      	inner_n = lengths_p[run_index] * !isna;
      	effective_width += inner_n;
      	temp_sum += values_p[run_index] * inner_n;
      }
      // Last run
      isna = ISNA(values_p[upper_run]);
      inner_n = (2 + run_first_index[upper_run] - (start + width)) * !isna;  //  +2 is because should also have ((start + width) - 1)
      effective_width += inner_n;
      temp_sum += values_p[upper_run] * inner_n;
      if ( effective_width != width && (effective_width == 0 || keep_na)) {
	temp_sum = na_val;
      } else {
	temp_sum /= effective_width;
      }
      ans_p[i] = temp_sum;
    }
  }
  UNPROTECT(1);
  return Ans;
}
   
  
SEXP RleViews_viewMeans3(SEXP Start, SEXP Width, SEXP Values, SEXP Lengths, SEXP Na_rm) {
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
  
  double temp_sum;
  int i, start, width, inner_n, effective_width;
  int lower_run, upper_run, run_index, mflag = 0;
  
  double* run_first_index = (double *) R_alloc(nrun, sizeof(double));
  double* run_last_index = (double *) R_alloc(nrun, sizeof(double));
  widthToStartEnd(lengths_p, run_first_index, run_last_index, nrun);
  
  // Abstract all the NA checking to a simple lookup of a boolean value
  char* isna = (char *) R_alloc(nrun, sizeof(char));
  isNA(Values, isna);

  // From here down all type-dependence could be handled by a template on values_p and na_val
  for (i = 0; i < nranges; i++) {
    start = start_p[i];
    width = width_p[i];
    temp_sum = 0;
    effective_width = 0;
    // Find run(s) covered by current range using something like findOverlaps(IRanges(start,width), ranges(rle))
    lower_run = findInterval(run_first_index, nrun, start, 0, 0, lower_run, &mflag) - 1;
    upper_run = findInterval(run_first_index, nrun, (start + width) - 1, 0, 0, lower_run, &mflag) - 1;  // Yes, search the left bound both times
    if (lower_run == upper_run) {  // Range all in one run, special case here allows simpler logic below
      ans_p[i] = values_p[lower_run];
      continue;
    } else {
      // First run
      inner_n = (1 + (run_last_index[lower_run] - start)) * !isna[lower_run];
      effective_width += inner_n;
      temp_sum += values_p[lower_run] * inner_n;
      // Inner runs
      for (run_index = lower_run + 1; run_index < upper_run; run_index++) {
      	inner_n = lengths_p[run_index] * !isna[run_index];
      	effective_width += inner_n;
      	temp_sum += values_p[run_index] * inner_n;
      }
      // Last run
      inner_n = (2 + run_first_index[upper_run] - (start + width)) * !isna[upper_run];  //  +2 is because should also have ((start + width) - 1)
      effective_width += inner_n;
      temp_sum += values_p[upper_run] * inner_n;
      if ( effective_width != width && (effective_width == 0 || keep_na)) {
	temp_sum = na_val;
      } else {
	temp_sum /= effective_width;
      }
      ans_p[i] = temp_sum;
    }
  }
  UNPROTECT(1);
  return Ans;
}
   
  
