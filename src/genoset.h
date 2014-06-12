#ifndef GENOSET_H
#define GENOSET_H

// utils.c
void isNA(SEXP vec, char* na);
int numNA(SEXP vec, char* na); // Would rather overload (C++ only) with int isNA(SEXP vec, char* na);
void widthToStartEnd(int* width, size_t* start, size_t* end, size_t n);
void widthToStart(int* width, size_t* start, size_t n);
size_t leftBound(size_t* values, size_t low, size_t high, size_t query);

#define LEFT_BOUND(values, low, high, query) do { \
  size_t probe = low + 1; \
  size_t jump = 2; \
  while (probe <= high && values[probe] <= query) { \
    low = probe; \
    probe += jump; \
    jump = jump << 1; \
  } \
  high = probe > high ? high + 1 : probe; \
  probe = ((high-low) >> 1) + low;   \
  while (low < probe) { \
    if (values[probe] > query) { \
      high = probe; \
    } else { \
      low = probe; \
    } \
    probe = ((high-low) >> 1) + low; \
    } \
  } while (0)
#endif

