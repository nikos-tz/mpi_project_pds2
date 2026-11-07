#include <stdio.h>
#include <stdlib.h>

double **alloc_2d_double(int rows, int cols){
  double *data = (double *)malloc(rows*cols*sizeof(double));
  double **array= (double **)malloc(rows*sizeof(double*));
  for(int i=0; i<rows; i++)
    array[i] = &(data[cols*i]);

  return array;
}

float **alloc_2d_float(int rows, int cols){
  float *data = (float *)malloc(rows*cols*sizeof(float));
  float **array= (float **)malloc(rows*sizeof(float*));
  for(int i=0; i<rows; i++)
    array[i] = &(data[cols*i]);

  return array;
}

int **alloc_2d_int(int rows, int cols){
  int *data = (int *)malloc(rows*cols*sizeof(int));
  int **array= (int **)malloc(rows*sizeof(int*));
  for(int i=0; i<rows; i++)
    array[i] = &(data[cols*i]);

  return array;
}