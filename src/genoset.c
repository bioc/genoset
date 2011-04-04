void binary_bound( int *starts, int *stops, int *data, int *initial_bounds, int *num_queries, int *start_indices, int *stop_indices ) {
  int query_index, probe, left, right, low, high;
  low = initial_bounds[0] - 2; /* Convert from 1-based to 0-based indices and set left of first allowable position, high is already there */

  for (query_index=0; query_index < num_queries[0]; query_index++) {
    high = initial_bounds[1];
    
    /* Left bound */
    left = starts[query_index];
    while (high - low > 1) {
      probe = (high + low) >> 1;
      if (data[probe] > left) {
        high = probe;
      } else {
	low = probe;
      }
    }
    start_indices[query_index] = low;
    
    /* Right bound */
    right = stops[query_index];
    high = initial_bounds[1];  /* Reset high to right end */
    while (high - low > 1) {
      probe = (high + low) >> 1;
      if (data[probe] < right) {
        low = probe;
      } else {
        high = probe;
      }
    }
    stop_indices[query_index] = high;
    low = start_indices[query_index]; /* Reset low to left end of this query to start next query */
    
  } /* End foreach query */
  
}
