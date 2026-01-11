#ifndef DISTRIBUTE_BY_MEDIAN
#define DISTRIBUTE_BY_MEDIAN

void distribute_by_median(int n, int d, int points_per_process,
                          int leader, int self_id, int num_tasks, int tag, MPI_Status status,
                          double** proc_data, double* pivot, double* distance_from_pivot);

#endif
