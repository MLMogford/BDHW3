#include <cstdio>
#include <math.h>
#include <cstdlib>
#include <sys/time.h>
#include <cstring>
#include <vector>
#include <thread>
#include <iostream>
//#include <mpif-sizeof.h>
#include<random>
#include <array>
#include <algorithm>
#include "mpi.h"
using namespace std;

//#include <boost/mpi.hpp>
#define REDUCE_REQUEST 1


#define MAX_LEN 4194304

long get_time_us()  // return the time in the unit of us
{
    struct timeval my_time;  //us
    gettimeofday(&my_time, nullptr);
    long runtime_us = 1000000 * my_time.tv_sec + my_time.tv_usec; // us
    return runtime_us;
}


/*
inline void test(size_t start, size_t end, size_t *result, int* a)//summation, passed as a parameter
{
    size_t sum = 0;
    for (size_t i=start; i<end; i++)
        sum += a[i];

    *result = sum;
}

void SUM_Threaded() {//thread launching
    std::vector<std::thread> threads(std::thread::hardware_concurrency());
    std::vector<size_t> results(std::thread::hardware_concurrency());

    size_t sum = 0;

    for (size_t j = 0; j < threads.size(); ++j) {
        size_t start =  j * a / std::thread::hardware_concurrency();
        size_t end = (j + 1) * a / std::thread::hardware_concurrency();
        threads[j] = std::thread(test, start, end, &results[j]);
    }

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < threads.size(); ++i) {
        threads[i].join();
        sum += results[i];
    }

    auto end = std::chrono::high_resolution_clock::now();

    cout << "Sum is " << sum << endl;

    int64_t elapse_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    cout<< "scan time: " << elapse_time << "usec" << endl;
    cout<< "Bandwidth is " << sizeof(int) * (double)1.0 * (a / elapse_time) <<" MB/s" << endl;
}

*/
int SUM_Array(int* a, int count) {
    size_t sum = 0;

    for (size_t i = 0; i < count; i++)
        sum += a[i];

    return sum;
}
//unthreaded test


   // auto start = std::chrono::high_resolution_clock::now();
    //test(0, *a, &sum);
   // auto end = std::chrono::high_resolution_clock::now();

  //  cout << "Sum is " << sum << endl;

  //  int64_t elapse_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
  //  cout<< "scan time: " << elapse_time << "usec" << endl;
  //  cout<< "Bandwidth is " << sizeof(int) * (double)1.0 * *a / elapse_time <<" MB/s" << endl;
//}



int main(int argc, char *argv[])
{
    std::vector<int> toto;
    int size, rank;

    MPI_Init(&argc,&argv);

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
    a = (int*)malloc(MAX_LEN * sizeof(int));            //origin array
    res = (int*)malloc(MAX_LEN * sizeof(int));
    res2 = (int*)malloc(MAX_LEN * sizeof(int));
    memset(a, 0 , MAX_LEN * sizeof(*a));                //sets to 0
    memset(res, 0 , MAX_LEN * sizeof(*res));
    memset(res2, 0 , MAX_LEN * sizeof(*res2));

    // TODO
    // you can add some variable or some other things as you want if needed
    // TODO

    for(count=4; count<=MAX_LEN; count*=16) // length of array : [ 4  64  1024  16384  262144  4194304 ]
    {
        // the element of array is generated randomly
        for(i=0; i<MAX_LEN; i++)
        {
            a[i] = rand() % MAX_LEN;
        }

        // MPI_Reduce and then print the usetime, the result will be put in res2[]
        MPI_Barrier(MPI_COMM_WORLD);
        begin_time = get_time_us();
        MPI_Reduce(a, res2, count, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        end_time = get_time_us();
        use_time2 = end_time - begin_time;
        if(rank == 0)
            printf("%d int use_time : %ld us [MPI_Reduce]\n", count, use_time2);


        // YOUR_Reduce and then print the usetime, the result should be put in res[]
        MPI_Barrier(MPI_COMM_WORLD);
        begin_time = get_time_us();
        // TODO

        if (rank > 0) { //not main process, send data to main
            MPI_Send(a, count, MPI_INTEGER, 0, REDUCE_REQUEST, MPI_COMM_WORLD);
        } else { //receive data from other process
            MPI_Status status;

            memcpy(res, a, count * sizeof(int));
            int source = 1;
            while (source < size) {
                MPI_Recv(a, count, MPI_INTEGER, MPI_ANY_SOURCE, REDUCE_REQUEST, MPI_COMM_WORLD, &status);
                for(i = 0; i < count; i++) {
                    res[i] += a[i];
                }

                source++;
            }
        }

        // you should delete the next line, and do the reduction using your idea


        // TODO
        MPI_Barrier(MPI_COMM_WORLD);
        end_time = get_time_us();
        use_time = end_time - begin_time;
        if(rank == 0)
            printf("%d int use_time : %ld us [YOUR_Reduce]\n", count, use_time);


        // check the result of MPI_Reduce and YOUR_Reduce
        if(rank == 0)
        {
            int correctness = 1;
            for(i=0; i<count; i++)
            {
                if(res2[i] != res[i])
                {
                    correctness = 0;
                }
            }
            if(correctness == 0)
                printf("WRONG !!!\n");
            else
                printf("CORRECT !\n");
        }


        if(rank == 0)
        {
            size_t sum = 0;
            begin_time = get_time_us();
            // TODO
            SUM_Array(a, count);
            //test - this is the function
            // TODO

  //          test(0, count, &sum);
    //        end_time = get_time_us();
      //      use_time = end_time - begin_time;
        //    printf("sum is %ld, use_time : %ld us [single thread]\n", sum, use_time);

            sum = 0;
            begin_time = get_time_us();
            // TODO
  //          SUM_Threaded();
            // calculate the sum of the result array reduced to process 0,
            // please make it faster with multiple threads.
            // TODO
            end_time = get_time_us();
            use_time = end_time - begin_time;
            printf("sum is %ld, use_time : %ld us [multiple threads]\n\n", sum, use_time);
        }

    }

    MPI_Finalize();

    return 0;
}
