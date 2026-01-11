#ifndef FUNCTIONS
#define FUNCTIONS

double **alloc_2d_double(int rows, int cols);

float **alloc_2d_float(int rows, int cols);

int **alloc_2d_int(int rows, int cols);

void merge(double* arr, int l, int m, int r);

void mergeSort(double* arr, int l, int r);

double find_max_double(double* arr, int length);

double find_min_double(double* arr, int length);

int trivial_log2_int (int x);

int trivial_pow_int(int base, int power);

#endif
