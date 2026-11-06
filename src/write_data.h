#ifndef WRITE_DATA_H
#define WRITE_DATA_H

#include <stdint.h>

/* Write the local points of an MPI process */
int write_proc_data(const char *directory, const char *filename, int n, int d, double **proc_data);

/* Write the pivot point */
int write_pivot(const char *filename, int d, double *pivot);

#endif