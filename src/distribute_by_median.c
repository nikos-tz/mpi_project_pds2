/**
 * Author: Nikos Tzatsis
 *
 * This function redistributes the points so that the first num_tasks/2 tasks
 * have the points with shorter or equal distance (from the pivot) to the median
 * distance and the remaining tasks have the longer or equal distances.
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "functions.h"
#include <sys/time.h>



void distribute_by_median(int n, int d, int points_per_process,
                          int leader, int self_id, int num_tasks, int tag, MPI_Status status,
                          double** proc_data, double* pivot, double* distance_from_pivot) {


    /**********/
    /* Get and send distance from pivot*/
    /**********/

    /**
     * Each task calculates the distances from the pivot and then sends them to the leader.
     * Note that in fact we calculate the square of the distance but this will not affect our output.
     */

    for(int i=0; i < points_per_process; ++i)
        distance_from_pivot[i] = 0.0;

    for(int i=0; i < points_per_process; ++i)
        for(int j=0; j < d; ++j)
            distance_from_pivot[i] += (pivot[j] - proc_data[i][j])*(pivot[j] - proc_data[i][j]);


    double myTime = 0.0;

    struct timeval start,end;


    double* all_distance_from_pivot = NULL;

    if (self_id != leader){
        MPI_Send(&(distance_from_pivot[0]), points_per_process, MPI_DOUBLE, leader, tag, MPI_COMM_WORLD);
    }
    else {

        all_distance_from_pivot = (double*) malloc(num_tasks * points_per_process * sizeof(double));

        for(int i=0; i < points_per_process; ++i)
            all_distance_from_pivot[i] = distance_from_pivot[i];

        double* sending_data;
        sending_data = (double*) malloc(points_per_process * sizeof(double));

        for(int i = leader + 1; i < (leader+num_tasks); ++i) {
            MPI_Recv(&(sending_data[0]), points_per_process, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);

            for(int j=0; j < points_per_process; ++j)
                all_distance_from_pivot[(i-leader)*points_per_process + j] = sending_data[j];
        }



        free(sending_data);
    }




    /**********/
    /* Calculate and send median distance */
    /**********/

    double median_distance;

    if(self_id == leader) {

        mergeSort(all_distance_from_pivot, 0, (num_tasks * points_per_process)-1);

        median_distance = (all_distance_from_pivot[(num_tasks * points_per_process)/2] + all_distance_from_pivot[((num_tasks * points_per_process)/2)-1])/2.0;


        for(int i = leader+1; i < (leader+num_tasks); ++i)
            MPI_Send(&median_distance, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
    }
    else {
        MPI_Recv(&median_distance, 1, MPI_DOUBLE, leader, tag, MPI_COMM_WORLD, &status);
    }






    /**********/
    /* order the distances*/
    /**********/

    /**
     * Each tasks orders its points internally. First go the points with the shorter than
     * the median distance and then the points with the equal one, and last the points with the longer one.
     * To do that we firstly give a value to every point depending its distance. Zero is for shorter,
     * one is for equal and two is for longer. We create a comparison array with this values.
     */
    int* comparison_arr;
    comparison_arr = (int*) calloc(points_per_process, sizeof(int));

    for(int i=0; i < points_per_process; ++i)
        if(distance_from_pivot[i] > median_distance)
            comparison_arr[i] = 2;
        else if(distance_from_pivot[i] == median_distance)
            comparison_arr[i] = 1;

    /**
     * Then we calculate the order of the points using the comparison array
     */

    int* ordered_indexes;
    ordered_indexes = (int*) malloc(points_per_process * sizeof(int));

    /**
     * This array contains the total number of the shorter, equal and longer
     * distances for this task
     */
    int* num_low_eq_high;
    num_low_eq_high = (int*) calloc(3, sizeof(int));

    int index = 0;

    for(int i=0; i < 3; ++i) {
        for(int j=0; j < points_per_process; ++j){

            if(comparison_arr[j] == i) {
                ordered_indexes[index] = j;
                ++index;
                ++num_low_eq_high[i];
            }

        }
    }

    /**
     * We use the prefix sum of the num_low_eq_high to reorder the points of each
     * task internally.
     */

    int* prefix_sum_own_low_eq_high;
    prefix_sum_own_low_eq_high = (int*) calloc(3, sizeof(int));

    for(int i=1; i < 3; ++i)
        prefix_sum_own_low_eq_high[i] = prefix_sum_own_low_eq_high[i-1] + num_low_eq_high[i-1];


    double** temp_proc_data;
    temp_proc_data = alloc_2d_double(points_per_process, d);

    for(int i=0; i < points_per_process; ++i)
        for(int j=0; j < d; ++j)
            temp_proc_data[i][j] = proc_data[ordered_indexes[i]][j];

    for(int i=0; i < points_per_process; ++i)
        for(int j=0; j < d; ++j)
            proc_data[i][j] = temp_proc_data[i][j];

    free(temp_proc_data[0]);
    free(temp_proc_data);




    /**********/
    /* make distribution matrices */
    /**********/

    /**
     * Each tasks sends its nom_low_eq_high array to the leader and then leader
     * connects them to create the all_low_eq_high matrix. The value all_low_eq_high[i][j]
     * tells us how many points with shorter (for i=0) or equal (for i=1) or longer (for i=2)
     * distance all the tasks before the task j have sent. If we add to this value the number
     * of all the points with shorter distance (for i=1) or the number of all the points with
     * shorter and equal distance (for i=2, we don't add something for i=0) and then we divide
     * the result with the number of the points per process (integer division) then we get the
     * receiver, that is, the task to which the task j must send its points which have shorter/
     * equal/longer distance (depending on i). If we don't divide that previews result but we
     * find the modular of the division above, then we find the place in the receiver, into which
     * the points must go. We accomplish the above complicate procedure in the next lines of code
     * by creating A LOT of prefix sum arrays.
     */

    int** all_low_eq_high;
    all_low_eq_high = alloc_2d_int(3, num_tasks);

    if(self_id == leader){
        for(int i=0; i < 3; ++i)
            all_low_eq_high[i][0] = num_low_eq_high[i];
    }

    for(int i=0; i < 3; ++i) {

        if(self_id != leader) {
            MPI_Send(&(num_low_eq_high[i]), 1, MPI_INT, leader, tag, MPI_COMM_WORLD);
        }
        else {

            for(int j = leader + 1; j < (leader+num_tasks); ++j) {
                MPI_Recv(&(all_low_eq_high[i][j-leader]), 1, MPI_INT, j, tag, MPI_COMM_WORLD, &status);
            }
        }

    }

    if(self_id == leader)
        for(int j = leader + 1; j < (leader+num_tasks); ++j)
            MPI_Send(&(all_low_eq_high[0][0]), 3*num_tasks, MPI_INT, j, tag, MPI_COMM_WORLD);
    else
        MPI_Recv(&(all_low_eq_high[0][0]), 3*num_tasks, MPI_INT, leader, tag, MPI_COMM_WORLD, &status);



    int* sum_all_low_eq_high;
    sum_all_low_eq_high = (int*) calloc(3, sizeof(int));

    int** all_prefix_sum_low_eq_high;
    all_prefix_sum_low_eq_high = alloc_2d_int(3, num_tasks);



    if(self_id == leader) {


        for(int i=0; i < 3; ++i)
            all_prefix_sum_low_eq_high[i][0] = 0;

        for(int i=0; i < 3; ++i)
            for(int j=1; j < num_tasks; ++j)
                all_prefix_sum_low_eq_high[i][j] = all_prefix_sum_low_eq_high[i][j-1] + all_low_eq_high[i][j-1];

        for(int i=0; i < 3; ++i)
            for(int j=0; j < num_tasks; ++j)
                sum_all_low_eq_high[i] += all_low_eq_high[i][j];


        for(int i = leader + 1; i < (leader+num_tasks); ++i)
                MPI_Send(&(all_prefix_sum_low_eq_high[0][0]), 3*num_tasks, MPI_INT, i, tag, MPI_COMM_WORLD);

        for(int i = leader + 1; i < (leader+num_tasks); ++i)
            MPI_Send(&(sum_all_low_eq_high[0]), 3, MPI_INT, i, tag, MPI_COMM_WORLD);



    }
    else {

        MPI_Recv(&(all_prefix_sum_low_eq_high[0][0]), 3*num_tasks, MPI_INT, leader, tag, MPI_COMM_WORLD, &status);

        MPI_Recv(&(sum_all_low_eq_high[0]), 3, MPI_INT, leader, tag, MPI_COMM_WORLD, &status);

    }


    int* prefix_sum_all_low_eq_high;
    prefix_sum_all_low_eq_high = (int*) calloc(3, sizeof(int));

    for(int i=1; i < 3; ++i)
        prefix_sum_all_low_eq_high[i] = prefix_sum_all_low_eq_high[i-1] + sum_all_low_eq_high[i-1];


    /**********/
    /* distribute the points */
    /**********/

    MPI_Request mpireq;



    temp_proc_data = alloc_2d_double(points_per_process, d);

    int wait_counter = 0;

    /**
     * Each task "watches" this for loop. If it has to send or receive points on some iteration then
     * it does it. If not, then it just skips the iteration. Note that if it has to send them to itself
     * then it doesn't send them via MPI but it simple copies them. Also we look up for the possibility
     * that one task will have to send its points to two continuously tasks (the receiver and the task
     * after it) because they dont fit in the receiver.
     */
    for(int i=0; i < 3; ++i) {


        for(int j=leader; j < (leader+num_tasks); ++j) {

            if(all_low_eq_high[i][j-leader] == 0) continue;

            int receiver = ( (prefix_sum_all_low_eq_high[i] + all_prefix_sum_low_eq_high[i][j-leader]) / points_per_process ) + leader;
            int position_in_receiver = (prefix_sum_all_low_eq_high[i] + all_prefix_sum_low_eq_high[i][j-leader]) % points_per_process;




            if((j+1-leader) < (num_tasks)) {

                int receiver_of_next = ( (prefix_sum_all_low_eq_high[i] + all_prefix_sum_low_eq_high[i][j+1-leader]) / points_per_process ) + leader;
                int position_in_receiver_of_next = (prefix_sum_all_low_eq_high[i] + all_prefix_sum_low_eq_high[i][j+1-leader]) % points_per_process;



                if( !( (receiver_of_next > receiver)&&(position_in_receiver_of_next > 0) ) ) {

                    if(receiver == j){
                        if(self_id == j){

                            int temp_index = position_in_receiver;
                            for(int k=prefix_sum_own_low_eq_high[i]; k < prefix_sum_own_low_eq_high[i]+all_low_eq_high[i][j-leader]; ++k){
                                for(int l=0; l < d; ++l){
                                        temp_proc_data[temp_index][l] = proc_data[k][l];
                                }
                                ++temp_index;
                            }
                        }
                    }
                    else if(self_id == j){
                        MPI_Send(&(proc_data[prefix_sum_own_low_eq_high[i]][0]), all_low_eq_high[i][j-leader]*d, MPI_DOUBLE, receiver, tag, MPI_COMM_WORLD);
                    }

                    else if(self_id == receiver){
                        MPI_Recv(&(temp_proc_data[position_in_receiver][0]), all_low_eq_high[i][j-leader]*d, MPI_DOUBLE, j, tag, MPI_COMM_WORLD, &status);
                    }
                }
                else {
                    int first_part = points_per_process - position_in_receiver;

                    if(receiver == j){
                        if(self_id == j){
                            int temp_index = position_in_receiver;
                            for(int k=prefix_sum_own_low_eq_high[i]; k < prefix_sum_own_low_eq_high[i]+first_part; ++k){
                                for(int l=0; l < d; ++l){
                                        temp_proc_data[temp_index][l] = proc_data[k][l];
                                }
                                ++temp_index;
                            }
                        }
                    }
                    else if(self_id == j){
                        MPI_Send(&(proc_data[prefix_sum_own_low_eq_high[i]][0]), first_part*d, MPI_DOUBLE, receiver, tag, MPI_COMM_WORLD);
                    }else if(self_id == receiver){
                        MPI_Recv(&(temp_proc_data[position_in_receiver][0]), first_part*d, MPI_DOUBLE, j, tag, MPI_COMM_WORLD, &status);
                    }
                    int second_part = all_low_eq_high[i][j-leader] - first_part;


                    if((receiver+1) == j){
                        if(self_id == j){
                            int temp_index = 0;
                            for(int k=prefix_sum_own_low_eq_high[i]+first_part; k < prefix_sum_own_low_eq_high[i]+first_part+second_part; ++k){
                                for(int l=0; l < d; ++l){
                                        temp_proc_data[temp_index][l] = proc_data[k][l];
                                }
                                ++temp_index;
                            }
                        }
                    }
                    else if(self_id == j){
                        MPI_Send(&(proc_data[prefix_sum_own_low_eq_high[i]+first_part][0]), second_part*d, MPI_DOUBLE, receiver+1, tag, MPI_COMM_WORLD);
                    }
                    else if(self_id == (receiver+1)){
                        MPI_Recv(&(temp_proc_data[0][0]), second_part*d, MPI_DOUBLE, j, tag, MPI_COMM_WORLD, &status);
                    }


                }

            }
            else {

                if(receiver == j){
                    if(self_id == j){

                        int temp_index = position_in_receiver;
                        for(int k=prefix_sum_own_low_eq_high[i]; k < prefix_sum_own_low_eq_high[i]+all_low_eq_high[i][j-leader]; ++k){
                            for(int l=0; l < d; ++l){
                                    temp_proc_data[temp_index][l] = proc_data[k][l];
                            }
                            ++temp_index;
                        }
                    }
                }
                else if(self_id == j){

                    MPI_Send(&(proc_data[prefix_sum_own_low_eq_high[i]][0]), all_low_eq_high[i][j-leader]*d, MPI_DOUBLE, receiver, tag, MPI_COMM_WORLD);


                }
                else if(self_id == receiver){
                    MPI_Recv(&(temp_proc_data[position_in_receiver][0]), all_low_eq_high[i][j-leader]*d, MPI_DOUBLE, j, tag, MPI_COMM_WORLD, &status);

                }
            }


       }

    }





    for(int i=0; i < points_per_process; ++i)
        for(int j=0; j < d; ++j)
            proc_data[i][j] = temp_proc_data[i][j];


    /**
     * Calculate the new distances from the pivot
     */

    for(int i=0; i < points_per_process; ++i)
        distance_from_pivot[i] = 0.0;

    for(int i=0; i < points_per_process; ++i)
        for(int j=0; j < d; ++j)
            distance_from_pivot[i] += (pivot[j] - proc_data[i][j])*(pivot[j] - proc_data[i][j]);



    /**********/
    /* free the pointers */
    /**********/

    if(self_id == leader)
        free(all_distance_from_pivot);

    free(comparison_arr);
    free(ordered_indexes);
    free(num_low_eq_high);
    free(all_low_eq_high[0]);
    free(all_low_eq_high);
    free(sum_all_low_eq_high);
    free(prefix_sum_all_low_eq_high);
    free(prefix_sum_own_low_eq_high);
    free(all_prefix_sum_low_eq_high[0]);
    free(all_prefix_sum_low_eq_high);
    free(temp_proc_data[0]);
    free(temp_proc_data);



}
