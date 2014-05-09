#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <R_ext/Utils.h>

// Take an SEXP and a pointer to a char array, fill out the char array 
//   with 0 (not NA) or 1 (NA) by looking at the values in vec
// Then we can compute on it without branching or polymorphism to account 
//   for type, Julia DataArrays style.
// This burns a bit of memory, but hopefully avoiding branching code for
//   NA checking will save a lot of time.
// Use the resulting NAs in a mathematical context rather than branching:
//   y += na[i] * x;
//   y += (! na[i]) * x;
// The compiler will pipeline the math for you, branches would prevent that.
// Argh, just noticed that Hadley thought of doing an isna C function, but his 
//   is inside out and returns INTSXP.
void isNA(SEXP vec, char* na) {
  if (TYPEOF(vec) == REALSXP) {
    double* vec_p = REAL(vec);    
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = ISNA(vec_p[i]);
    }
  } else if (TYPEOF(vec) == INTSXP) {
    int* vec_p = INTEGER(vec);
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = vec_p[i] == NA_INTEGER;
    }
  } else if (TYPEOF(vec) == LGLSXP) {
    int* vec_p = INTEGER(vec);
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = vec_p[i] == NA_LOGICAL;
    }
  } else if (TYPEOF(vec) == STRSXP) {
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = STRING_ELT(vec, i) == NA_STRING; 
    }
  } else {
    error("vec must contain either 'integer', 'logical' or 'character' or 'numeric' values");
  }
}
