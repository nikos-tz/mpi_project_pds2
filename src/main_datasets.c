/**
 * Author: Nikos Tzatsis
 * Title: Distribute points by median distance
 *
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "distribute_by_median.h"
#include "functions.h"



int main( int argc, char** argv ) {

    FILE* fp = NULL;
    /* insert here the path to your dataset file */
    fp = fopen("/home/csal/pds/pds-codebase/datasets/cifar10.bin", "r");

    /* number of points and dimension of points */
    int n, d;

    /* variables for the communicator */
    int tag = 1;
    int num_tasks;
    int self_id;
    int leader = 0; // leader task

    /* initialise the communicator */
    MPI_Status status;
    MPI_Init (&argc,&argv);
    MPI_Comm_size (MPI_COMM_WORLD,&num_tasks);
    MPI_Comm_rank (MPI_COMM_WORLD,&self_id);

    int64_t n_file, d_file;

    int ignore = 0; //just to ignore the return values of some functions

    /**
     * The leader reads the dataset file. At first, it gets the N,D from the file, displays them
     * and then lets the user to enter their own n,d values which must be powers of 2
     * and also less than N,D. Then send n,d to the other tasks.
     */
    if(self_id == leader){

        ignore = fread(&d_file, sizeof(int64_t), 1, fp);
        ignore = fread(&n_file, sizeof(int64_t), 1, fp);

        printf("N = %ld and D = %ld\n", n_file, d_file);

        printf("Enter n: \n");

        ignore = scanf("%d", &n);

        printf("\nEnter d: \n");

        ignore = scanf("%d", &d);

        for(int i=leader+1; i < (leader+num_tasks); ++i)
            MPI_Send(&n, 1, MPI_INT, i, tag, MPI_COMM_WORLD);

        for(int i=leader+1; i < (leader+num_tasks); ++i)
            MPI_Send(&d, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
    }
    else {
        MPI_Recv(&n, 1, MPI_INT, leader, tag, MPI_COMM_WORLD, &status);

        MPI_Recv(&d, 1, MPI_INT, leader, tag, MPI_COMM_WORLD, &status);
    }

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


        for(int i=0; i < points_per_process; ++i) {
            /* we ignore some values from the file, because our n or d could be less than the N,D of the file */
            fseek(fp, (2*sizeof(int64_t) + (leader*points_per_process*d_file + i*d_file)*sizeof(double)), SEEK_SET);

            ignore = fread(&proc_data[i][0], sizeof(double), d, fp);

        }

        for(int i=leader+1; i < (leader+num_tasks); ++i) {
            for(int j=0; j < points_per_process; ++j) {

                fseek(fp, (2*sizeof(int64_t) + (i*points_per_process*d_file + j*d_file)*sizeof(double)), SEEK_SET);

                ignore = fread(&sending_data[j][0], sizeof(double), d, fp);
            }

            MPI_Send(&(sending_data[0][0]), d*points_per_process, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);

        }

    }
    else{
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

    }
    else {
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


    for(int i=0; i < trivial_log2_int(num_tasks); ++i) {

        for(int j=0; j < trivial_pow_int(2, i); ++j) {

            int current_num_tasks = num_tasks / trivial_pow_int(2, i);
            int current_leader = j * current_num_tasks;

            if( (self_id >= current_leader)&&(self_id < (current_leader + current_num_tasks)) ) {
                distribute_by_median(n, d, points_per_process,
                                        current_leader, self_id, current_num_tasks, tag, status,
                                        proc_data, pivot, distance_from_pivot);
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

    if(self_id == leader) {
        all_max_distance = (double*) malloc(num_tasks*sizeof(double));
        all_min_distance = (double*) malloc(num_tasks*sizeof(double));
    }

    if(self_id != leader) {
        MPI_Send(&max_distance, 1, MPI_DOUBLE, leader, tag, MPI_COMM_WORLD);
    }
    else{
        all_max_distance[0] = max_distance;

        for(int i=leader+1; i < (leader+num_tasks); ++i){
            MPI_Recv(&all_max_distance[i-leader], 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
        }
    }

    if(self_id != leader) {
        MPI_Send(&min_distance, 1, MPI_DOUBLE, leader, tag, MPI_COMM_WORLD);
    }
    else {
        all_min_distance[0] = min_distance;

        for(int i=leader+1; i < (leader+num_tasks); ++i){
            MPI_Recv(&all_min_distance[i-leader], 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
        }

        int error_counter = 0;

        for(int i=1; i < num_tasks; ++i){
            if(all_min_distance[i] < all_max_distance[i-1])
                ++error_counter;
        }

        printf("---------------DISTANCES---------------\n");
        printf("TASK\t\tMIN\t\tMAX\n\n");

        for(int i=0; i < num_tasks; ++i)
            printf("%d\t\t%lf\t%lf\n", i, all_min_distance[i], all_max_distance[i]);


        printf("\nNumber of errors is: %d\n\n", error_counter);

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


