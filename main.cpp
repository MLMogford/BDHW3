#include <cstdio>
#include <math.h>
#include <cstdlib>
#include <sys/time.h>
#include <cstring>
#include <vector>
#include <thread>
#include <iostream>
#include<random>
#include <array>
#include <algorithm>
#include "mpi.h"

using namespace std;

#define REDUCE_REQUEST 1


#define MAX_LEN 4194304

long get_time_us()  // return the time in the unit of us
{
    struct timeval my_time;  //us
    gettimeofday(&my_time, nullptr);
    long runtime_us = 1000000 * my_time.tv_sec + my_time.tv_usec; // us
    return runtime_us;
}



int SUM_Array(int *a, int count) {
    int sum = 0;

    for (size_t i = 0; i < count; i++)
        sum += a[i];

    return sum;
}
//unthreaded test

int SUM_Array_threaded(int *a, int count) {
    int sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (size_t i = 0; i < count; i++)
        sum += a[i];

    return sum;
}





int main(int argc, char *argv[]) {
    std::vector<int> toto;
    int size, rank;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(time(NULL));  // seed for rand() to generate the array randomly

    int *a;     // the array used to do the reduction
    int *res;   // the array to record the result of YOUR_Reduce   - data is in res
    int *res2;  // the array to record the result of MPI_Reduce

    int count;
    long begin_time, end_time, use_time, use_time2; // use_time for YOUR_Reduce & use_time2 for MPI_Reduce

    int i;

    // initialize
    a = (int *) malloc(MAX_LEN * sizeof(int));            //origin array
    res = (int *) malloc(MAX_LEN * sizeof(int));
    res2 = (int *) malloc(MAX_LEN * sizeof(int));
    memset(a, 0, MAX_LEN * sizeof(*a));                //sets to 0
    memset(res, 0, MAX_LEN * sizeof(*res));
    memset(res2, 0, MAX_LEN * sizeof(*res2));

    // TODO
    // you can add some variable or some other things as you want if needed
    // TODO

    volatile int force_create_thread = SUM_Array_threaded(a, MAX_LEN);

    for (count = 4; count <= MAX_LEN; count *= 16) // length of array : [ 4  64  1024  16384  262144  4194304 ]
    {
        // the element of array is generated randomly
        for (i = 0; i < MAX_LEN; i++) {
            a[i] = rand() % MAX_LEN;
        }

        // MPI_Reduce and then print the usetime, the result will be put in res2[]
        MPI_Barrier(MPI_COMM_WORLD);
        begin_time = get_time_us();
        MPI_Reduce(a, res2, count, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        end_time = get_time_us();
        use_time2 = end_time - begin_time;
        if (rank == 0)
            printf("%d int use_time : %ld us [MPI_Reduce]\n", count, use_time2);


        // YOUR_Reduce and then print the usetime, the result should be put in res[]
        MPI_Barrier(MPI_COMM_WORLD);
        begin_time = get_time_us();
        /*/ TODO

        if (rank > 0) { //not main process, send data to main
            MPI_Send(a, count, MPI_INTEGER, 0, REDUCE_REQUEST, MPI_COMM_WORLD);
        } else { //receive data from other process
            MPI_Status status;

            memcpy(res, a, count * sizeof(int));
            int source = 1;
            while (source < size) {
                MPI_Recv(a, count, MPI_INTEGER, MPI_ANY_SOURCE, REDUCE_REQUEST, MPI_COMM_WORLD, &status);
                for (i = 0; i < count; i++) {
                    res[i] += a[i];
                }

                source++;
            }
        }
*/


        MPI_Barrier(MPI_COMM_WORLD);
        begin_time = get_time_us();
        // TODO

        if (rank %2 != 0) { //odd no. process, send data to even
            MPI_Send(a, count, MPI_INTEGER, rank - 1, REDUCE_REQUEST, MPI_COMM_WORLD);

        } else if (rank %2 ==0 && rank != 4 && rank != 0) {
            MPI_Send(a, count, MPI_INTEGER, rank -2, REDUCE_REQUEST, MPI_COMM_WORLD);
        }
        else  if (rank  == 4) {
            MPI_Send(a, count, MPI_INTEGER, 0, REDUCE_REQUEST, MPI_COMM_WORLD);
        }

         else { //receive data from other process
            MPI_Status status;
            int result;
            result = log2(size);
            memcpy(res, a, count * sizeof(int));
            for (int j = 0 ; j < result ; j++) {
                MPI_Recv(a, count, MPI_INTEGER, MPI_ANY_SOURCE, REDUCE_REQUEST, MPI_COMM_WORLD, &status);
                for (i = 0; i < count; i++) {
                    res[i] += a[i];
                }
            }
        }


        // TODO
        MPI_Barrier(MPI_COMM_WORLD);
        end_time = get_time_us();
        use_time = end_time - begin_time;
        if (rank == 0)
            printf("%d int use_time : %ld us [YOUR_Reduce]\n", count, use_time);


        // check the result of MPI_Reduce and YOUR_Reduce
        if (rank == 0) {
            int correctness = 1;
            for (i = 0; i < count; i++) {
                if (res2[i] != res[i]) {
                    correctness = 0;
                }
            }
            if (correctness == 0)
                printf("WRONG !!!\n");
            else
                printf("CORRECT !\n");
        }


        if (rank == 0) {
            int sum = 0;
            begin_time = get_time_us();
            sum = SUM_Array(res, count);
            end_time = get_time_us();
            use_time = end_time - begin_time;
            printf("sum is %i, use_time : %ld us [single thread]\n", sum, use_time);

            sum = 0;
            begin_time = get_time_us();
            sum = SUM_Array_threaded(res, count);
            end_time = get_time_us();
            use_time = end_time - begin_time;
            printf("sum is %i, use_time : %ld us [multiple threads]\n\n", sum, use_time);
        }

    }

    MPI_Finalize();

    return 0;
}
