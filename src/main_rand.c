/**
 * Author: Nikos Tzatsis
 * Title: Distribute points by median distance
 *
 *
 * This function is exactly like the main_datasets function, but here
 * our data are randomly created.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "distribute_by_median.h"
#include "functions.h"



int main( int argc, char** argv ) {

    int n, d;

    int tag = 1;
    int num_tasks;
    int self_id;
    int leader = 0;

    MPI_Status status;
    MPI_Init (&argc,&argv);
    MPI_Comm_size (MPI_COMM_WORLD,&num_tasks);
    MPI_Comm_rank (MPI_COMM_WORLD,&self_id);

    int ignore = 0;

    if(self_id == leader){


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

    int points_per_process = n/num_tasks;

    /*********/
    /* for random data*/
    /*********/

    double **point = NULL;

    if(self_id == leader){
        point = alloc_2d_double(n,d);


        srand((unsigned int)time(NULL));

        double a = 10.0;

        for (int i=0; i < n; ++i)
            for (int j=0; j < d; ++j)
                point[i][j] = ( (double)rand() / (double)(RAND_MAX) ) * a;
    }

    double** proc_data;
    proc_data = alloc_2d_double(points_per_process, d);

    double** sending_data = NULL;


    /*********/
    /* for random data*/
    /*********/

    if (self_id == leader) {

        sending_data = alloc_2d_double(points_per_process, d);


        for ( int i=0; i < points_per_process; ++i)
            for (int j=0; j < d; ++j)
                proc_data[i][j] = point[i][j];



        for (int i=1; i < num_tasks; ++i) {

            for (int k=0; k < points_per_process; ++k)
                for (int j=0; j < d; ++j)
                    sending_data[k][j] = point[i*points_per_process + k][j];


            MPI_Send(&(sending_data[0][0]), d*points_per_process, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
        }

    }
    else {
        MPI_Recv(&(proc_data[0][0]), d*points_per_process, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }



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






    if(self_id == leader){
        free(point[0]);
        free(point);
    }



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


