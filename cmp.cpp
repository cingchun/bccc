/*****************************************************
      Copyright:   Qingchun Wang @ NJU
      File Name:   cmp.c
            Des:   c multiprocess(parallel)
           Mail:   qingchun720@foxmail.com
   Created Time:   16:48 10月-22/2019
******************************************************/

#include <iostream>
#include <ctime>
#include "stdio.h"
#include "omp.h"

using namespace std;

long long N=100000;

extern "C"
{
    long long cmp(void)
    {
        int nthread,ithread;
        long long n = 0.0;

        /* test parallel */
        //#pragma omp parallel private(nthread,ithread)
        //{
        //    //obtian and print thread id
        //    ithread=omp_get_thread_num();
        //    printf("Hello Word from OMP0 thread %d\n",ithread);
        //    // only master thread does this;
        //    if(ithread==0)
        //    {
        //        nthread = omp_get_num_threads();
        //        printf("Number of thread: %d\n", nthread);
        //    }
        //}
        /* test sections */
        //#pragma omp parallel sections private(nthread,ithread)
        //{
        //    #pragma omp section
        //    {
        //        //obtian and print thread id
        //        ithread=omp_get_thread_num();
        //        printf("Hello Word from OMP0 thread %d\n",ithread);
        //    }
        //    #pragma omp section
        //    {
        //        //obtian and print thread id
        //        ithread=omp_get_thread_num();
        //        printf("Hello Word from OMP1 thread %d\n",ithread);
        //    }
        //    #pragma omp section
        //    {
        //        //obtian and print thread id
        //        ithread=omp_get_thread_num();
        //        printf("Hello Word from OMP2 thread %d\n",ithread);
        //    }
        //}
        /* test for */
        #pragma omp parallel for //num_threads(1)
        for (long long i=0; i<N; ++i)
        {
            //ithread=omp_et_thread_num();
            //printf("from OMP thread %d: i = %d\n",ithread, i);
            long long sum = 0;
            for (long long j=0; j<N; ++j)
                for (long long k=0; k<N; ++k)
                    sum += (i+1)*(j+1)*(k+1);
            #pragma omp atomic
            n += sum;
        }

       return n;
    }
}



int main(void)
{
    long long n;
    time_t t1, t2;
	printf("*** Test python call parallel c ***\n");
	time(&t1);
	n = cmp();
	printf(" n = %ld \n", n);
	time(&t2);
	printf(" Elapsed time: %.6f sec. \n", difftime(t2,t1));
	printf("*** cloze testing ***\n");
	return 0;
}




