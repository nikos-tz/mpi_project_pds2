/**
 * Title: Distribute points by median distance
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "distribute_by_median.h"
#include <src/util/output.h>
#include <src/util/file.h>
#include <src/aux/memory.h>
#include <src/aux/math.h>
#include <src/aux/sort.h>

void check_input(const int argc, char** argv){
  if(argc != 2){
    fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
}

void initialise_communicator(int argc, char** argv, int* num_tasks, int* self_id){
  MPI_Init(&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, num_tasks);
  MPI_Comm_rank (MPI_COMM_WORLD, self_id);
}

void dataset_size(const int self_id, const int leader, FILE* fp, int64_t* n_file, int64_t* d_file, int* n, int* d){
  /**
   * The leader reads the dataset file. At first, it gets the N,D from the file, displays them
   * and then lets the user to enter their own n,d values which must be powers of 2
   * and also less than N,D. Then send n,d to the other tasks.
   */
  if(self_id == leader){
    file_read_int64(d_file, 1, fp);
    file_read_int64(n_file, 1, fp);

    printf("N = %ld and D = %ld\n", *n_file, *d_file);

    printf("Enter n: \n");

    scan_one_int(n);

    printf("\nEnter d: \n");

    scan_one_int(d);
  }
  MPI_Bcast(n, 1, MPI_INT, leader, MPI_COMM_WORLD);
  MPI_Bcast(d, 1, MPI_INT, leader, MPI_COMM_WORLD);
}

int main(int argc, char** argv){
  int ignore = 0; //just to ignore the return values of some functions

  check_input(argc, argv);

  FILE* fp = open_file(argv[1]);

  /* variables for the communicator */
  MPI_Status status;
  int tag = 1;
  int leader = 0;
  int num_tasks;
  int self_id;

  initialise_communicator(argc, argv, &num_tasks, &self_id);

  /* number of points and dimension of points */
  int64_t n_file, d_file;
  int n, d;
  dataset_size(self_id, leader, fp, &n_file, &d_file, &n, &d);

  int points_per_process = n/num_tasks; // how many points each process will hold

  /* declare a 2D array for our points, each row is a point */
  double** proc_data;
  proc_data = alloc_2d_double(points_per_process, d);

  double** sending_data = NULL;

  /**
   * Now leader reads the points from the file and then sends them to the other tasks
   */
  if(self_id == leader) {
    sending_data = alloc_2d_double(points_per_process, d);

    for(int i=0; i < points_per_process; ++i){
      /* we ignore some values from the file, because our n or d could be less than the N,D of the file */
      fseek(fp, (2*sizeof(int64_t) + (leader*points_per_process*d_file + i*d_file)*sizeof(double)), SEEK_SET);

      ignore = fread(&proc_data[i][0], sizeof(double), d, fp);
    }

    for(int i=leader+1; i < (leader+num_tasks); ++i){
      for(int j=0; j < points_per_process; ++j){
        fseek(fp, (2*sizeof(int64_t) + (i*points_per_process*d_file + j*d_file)*sizeof(double)), SEEK_SET);
        ignore = fread(&sending_data[j][0], sizeof(double), d, fp);
      }

      MPI_Send(&(sending_data[0][0]), d*points_per_process, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(&(proc_data[0][0]), d*points_per_process, MPI_DOUBLE, leader, tag, MPI_COMM_WORLD, &status);
  }

  fclose(fp);

  /*********/
  /* Get and send the pivot*/
  /*********/

  double* pivot = NULL;

  pivot = (double*) malloc(d * sizeof(double));

  if(self_id == leader){
    for(int i=0; i < d; ++i)
      pivot[i] = proc_data[0][i];

    for(int i = leader + 1; i < (leader+num_tasks); ++i)
      MPI_Send(&(pivot[0]), d, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
  } else {
    MPI_Recv(&(pivot[0]), d, MPI_DOUBLE, leader, tag, MPI_COMM_WORLD, &status);
  }

  /**
   * We call distribute_by_median() multiple times in order to distribute points
   * which are divided into subdivisions of distinct distances in each process
   * (that is, the first process will have the points with the shortest distances,
   * the second process will have the points with the second shortest distances etc).
   * To do that, the function calls are made in log2(p) levels (we will count them
   * with variable i), where p is the total number of tasks. In each level we call
   * the function 2^i times (we count those times with variable j). Each of this times,
   * the function is called with (p / 2^i) number of tasks and j * (p / 2^i) as the leader.
   */
  double* distance_from_pivot;
  distance_from_pivot = (double*) malloc(points_per_process * sizeof(double));

  for(int i=0; i < trivial_log2_int(num_tasks); ++i){
    for(int j=0; j < trivial_pow_int(2, i); ++j){
      int current_num_tasks = num_tasks / trivial_pow_int(2, i);
      int current_leader = j * current_num_tasks;

      if( (self_id >= current_leader)&&(self_id < (current_leader + current_num_tasks)) ){
        distribute_by_median(
          n, d, points_per_process, current_leader, self_id, current_num_tasks, tag, status, proc_data, pivot, distance_from_pivot
        );
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }

    if(self_id == leader){
      printf("\n");
    }
  }

  /**
   * Now we check if the correctness of our results. To do that, the leader collects
   * the min and max distance of each task. Then we check if the max distance of a
   * task is greater than the min distance of the next task. If that happens, an error
   * counter is increased. At the the end we display all the min and max distances and
   * the error counter, which must be equal to zero.
   */

  double max_distance = find_max_double(distance_from_pivot, points_per_process);
  double min_distance = find_min_double(distance_from_pivot, points_per_process);

  double* all_max_distance = NULL;
  double* all_min_distance = NULL;

  if(self_id == leader){
    all_max_distance = (double*) malloc(num_tasks*sizeof(double));
    all_min_distance = (double*) malloc(num_tasks*sizeof(double));
  }

  if(self_id != leader){
      MPI_Send(&max_distance, 1, MPI_DOUBLE, leader, tag, MPI_COMM_WORLD);
  } else {
    all_max_distance[0] = max_distance;

    for(int i=leader+1; i < (leader+num_tasks); ++i){
      MPI_Recv(&all_max_distance[i-leader], 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
    }
  }

  if(self_id != leader){
    MPI_Send(&min_distance, 1, MPI_DOUBLE, leader, tag, MPI_COMM_WORLD);
  } else {
    all_min_distance[0] = min_distance;

    for(int i=leader+1; i < (leader+num_tasks); ++i){
      MPI_Recv(&all_min_distance[i-leader], 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
    }

    int error_counter = 0;

    for(int i=1; i < num_tasks; ++i){
      if(all_min_distance[i] < all_max_distance[i-1]) ++error_counter;
    }

    printf("---------------DISTANCES---------------\n");
    printf("TASK\t\tMIN\t\tMAX\n\n");

    for(int i=0; i < num_tasks; ++i)
      printf("%d\t\t%lf\t%lf\n", i, all_min_distance[i], all_max_distance[i]);

    printf("\nNumber of errors is: %d\n\n", error_counter);
  }

  /* Write each process's local data */

  if(self_id == leader){
    if(delete_directory("proc_data") != 0){
      fprintf(stderr, "Could not delete proc_data directory\n");
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  char filename[100];

  sprintf(filename, "proc_data_rank_%d.bin", self_id);

  if(write_proc_data("proc_data", filename, points_per_process, d, proc_data) != 0){
    printf("Rank %d: Error writing proc_data\n", self_id);
  }

  /* Only the leader writes the pivot */

  if(self_id == leader){
    if(write_pivot("pivot.bin", d, pivot) != 0){
      printf("Error writing pivot\n");
    }
  }

  /* free all the dynamically allocated variables  */

  if(self_id == leader){
    free(sending_data[0]);
    free(sending_data);

    free(all_max_distance);
    free(all_min_distance);
  }

  free(proc_data[0]);
  free(proc_data);
  free(pivot);
  free(distance_from_pivot);

  MPI_Finalize();

  return (0);
}


