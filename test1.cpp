#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>

int main(int argc, char *argv[]){
  double time1=0.0, tstart;
  int nthreads, tid;
  tstart = clock();
  time1 += clock() - tstart;
  time1 = time1/CLOCKS_PER_SEC;
  printf("%lf\n",time1 );

  #pragma omp parallel for
  for (size_t i = 0; i < 1000; i++) {
    printf("I'm thread %i of %d\n", omp_get_thread_num(), omp_get_num_threads());
  }
  /* Fork a team of threads giving them their own copies of variables */
  #pragma omp parallel private(nthreads, tid)
    {
    /* Obtain thread number */
    tid = omp_get_thread_num();
    printf("Hello World from thread = %d\n", tid);

    /* Only master thread does this */
    if (tid == 0)
      {
      nthreads = omp_get_num_threads();
      printf("Number of threads = %d\n", nthreads);
      }



    }  /* All threads join master thread and disband */

}
