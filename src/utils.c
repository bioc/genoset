#include <R.h>
#include <Rinternals.h>
#include "genoset.h"

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
//  This could alternatively be implemented as vec<bool>. Probably very similar speed, but smaller.
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

// Like the above isNA, but returns the number of NAs.
//   This may be used to decide proceed with a simpler algorithm 
//   if == 0 or the sum itself may be used to simplify subsequent
//   code, like a mean.
int numNA(SEXP vec, char* na) {
  int num_na = 0;
  if (TYPEOF(vec) == REALSXP) {
    double* vec_p = REAL(vec);    
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = ISNA(vec_p[i]);
      num_na += na[i];
    }
  } else if (TYPEOF(vec) == INTSXP) {
    int* vec_p = INTEGER(vec);
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = vec_p[i] == NA_INTEGER;
      num_na += na[i];
    }
  } else if (TYPEOF(vec) == LGLSXP) {
    int* vec_p = INTEGER(vec);
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = vec_p[i] == NA_LOGICAL;
      num_na += na[i];
    }
  } else if (TYPEOF(vec) == STRSXP) {
    for (int i=0; i < LENGTH(vec); i++) {
      na[i] = STRING_ELT(vec, i) == NA_STRING; 
      num_na += na[i];
    }
  } else {
    error("vec must contain either 'integer', 'logical' or 'character' or 'numeric' values");
  }
  return(num_na);
}

// Take widths, like from an Rle, compute start and end index for each run (1-based).
// Want to switch to unsigned int as it will fit a human genome pos.
void widthToStartEnd(int* width, int* start, int* end, int n) {
  start[0] = 1;
  end[0] = width[0];
  for (int i=1; i < n; i++) {
    start[i] = end[i-1] + 1;
    end[i] = end[i-1] + width[i];
  }
}

// Take widths, like from an Rle, compute start index for each run (1-based).
void widthToStart(int* width, int* start, int n) {
  start[0] = 1;
  for (int i=1; i < n; i++) {
    start[i] = start[i-1] + width[i-1];
  }
}

/***********************
Searching
***********************/
// unsigned int is big enough for a position in the human genome, but 32 bits
// Also consider using unsigned ints to save having to care negative values passed in
// Fun fact: if you do --positions before you pass it in, you get 1-based indices back, a al findInterval.
int leftBound(int* positions, int query, int n, unsigned int restart) {
  unsigned int probe, low, high, gap, jump;
  // If data unsorted, current target may be anywhere left of restart for previous target, just start at 0
  //  printf("\nleftBound. low: %i, high: %i, positions[high]: %i, n: %i\n", low, high, positions[high], n);
  high = low = restart < n ? restart : n - 1;
  //  printf("leftBound. low: %i, high: %i, positions[high]: %i, n: %i\n", low, high, positions[high], n);
  if (query < positions[low]) {
    low = 0;  // Consider making this -1 for generality
  }
  //  printf("leftBound. low: %i, high: %i, positions[high]: %i, n: %i\n", low, high, positions[high], n);
  // Right bound likely close to previous
  for (jump=1;  ;jump+=jump) {
    high += jump;
    if (high >= n) {
      high = n;
      break;
    }
    if (query < positions[high]) {
      break;
    }
    low = high;
  }
  //  printf("leftBound. low: %i, high: %i, positions[high]: %i, n: %i\n", low, high, positions[high], n);
  // Now binary search for closest left bound
  gap = high - low;
  while (gap > 1) {
    probe = (gap >> 1) + low;
    if (positions[probe] > query) {
      high = probe;
    } else {
      low = probe;
    }
    gap = high - low;  // Hack to avoid integer overflow with low + high
  }
  //  printf("leftBound. low: %i, high: %i, positions[high]: %i, n: %i\n\n", low, high, positions[high], n);
  return(low);
}

//int leftBound2(int* first int* last, int query, int* probe) {
//  unsigned int gap, jump;
//  // If data unsorted, current target may be anywhere left of restart for previous target, just start at 0
//  //  printf("\nleftBound. low: %i, high: %i, positions[high]: %i, n: %i\n", low, high, positions[high], n);
//  probe = probe < last ? probe : last;
//  probe = probe >= first && *probe > *first ? probe : first;
//  //  printf("leftBound. low: %i, high: %i, positions[high]: %i, n: %i\n", low, high, positions[high], n);
//
//  // Right bound likely close to previous
//  for (jump=1;  ;jump+=jump) {
//    high += jump;
//    if (high >= n) {
//      high = n;
//      break;
//    }
//    if (query < positions[high]) {
//      break;
//    }
//    low = high;
//  }
//  //  printf("leftBound. low: %i, high: %i, positions[high]: %i, n: %i\n", low, high, positions[high], n);
//  // Now binary search for closest left bound
//  gap = high - low;
//  while (gap > 1) {
//    probe = (gap >> 1) + low;
//    if (positions[probe] > query) {
//      high = probe;
//    } else {
//      low = probe;
//    }
//    gap = high - low;  // Hack to avoid integer overflow with low + high
//  }
//  //  printf("leftBound. low: %i, high: %i, positions[high]: %i, n: %i\n\n", low, high, positions[high], n);
//  return(low);
//}
