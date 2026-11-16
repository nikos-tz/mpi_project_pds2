/**
 * Title: Distribute points by median distance
 *
 *
 * This function is exactly like the main_datasets function, but here
 * our data are randomly created.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <src/util/output.h>
#include <src/util/file.h>
#include <src/aux/memory.h>

int main(int argc, char** argv){
  int n, d;

  int ignore = 0;

  printf("Enter n: \n");

  ignore = scanf("%d", &n);

  printf("\nEnter d: \n");

  ignore = scanf("%d", &d);

  double **point = NULL;

  point = alloc_2d_double(n,d);

  srand((unsigned int)time(NULL));

  double a = 10.0;

  for (int i=0; i < n; ++i)
    for (int j=0; j < d; ++j)
      point[i][j] = ( (double)rand() / (double)(RAND_MAX) ) * a;

  if(write_proc_data("datasets", "rand.bin", n, d, point) != 0){
    printf("Error writing proc_data\n");
  }

  return (0);
}


