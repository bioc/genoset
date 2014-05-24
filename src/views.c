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
   int i, start, width, index, lower_run, upper_run, lower_bound, upper_bound, max_index;
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
     // I think doing genoset::binary_bound on the cumsums of start and end would be faster here for finding lower and upper run.
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
  
// Some branching for NA check only necessary for na.rm=TRUE and numeric Values.
// Could let the NaNs contaminate for na.rm=FALSE or use x * !na[i] for int types, 
//   but this is a minor fraction of the time. 75% of time in findInterval. Could 
//   make version that does not check target list for NAs (guaranteed
//   none as coming from rle runLength) gives 0-based result indices (hmm, --pointer 
//   optionally to switch?) and possibly using long ints for positions to search against.
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
