#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>

double function(double x, double y){
  return x * x * ( - (y * y));
}

double partialintegral(double hxmin, double hxmax, double hymin, double hymax, int nxmax, int nymax){
  double integral = 0;


double integrate(double hxmin, double hxmax, double hymin, double hymax, int nxmax, int nymax){
  double deltahx, deltahy;
  double integral = 0;
  double time1 = 0.0, tstart;
  double nstart;
  deltahx = (hxmax - hxmin) / nxmax;
  deltahy = (hymax - hymin) / nymax;
  tstart = clock();
  omp_get_num_threads()
  nystart = omp_get_thread_num() * nymax / omp_get_num_threads()
  for (int ny = nystart; ny < nymax; ny++) {
    double c = 0;
    for (int nx = 0; nx < nxmax; nx++) {
      c += deltahx * deltahy / 4 * (
        function(hxmin + deltahx * nx      , hymin + deltahy * ny      ) +
        function(hxmin + deltahx * (nx + 1), hymin + deltahy * ny      ) +
        function(hxmin + deltahx * nx      , hymin + deltahy * (ny + 1)) +
        function(hxmin + deltahx * (nx + 1), hymin + deltahy * (ny + 1)));
      //printf("c = %f\n", c);
      //printf("int = %f\n", integral);
  return integral

    }
    //integral +=c;
  }
  time1 += clock() - tstart;
  time1 = time1/CLOCKS_PER_SEC;
  printf("%f\n", time1);
  return integral;
}
int main(int argc, char *argv[]){
  double time1=0.0, tstart;
  int nthreads, tid;
  double hxmin = 0., hymin = 0., hxmax = 1., hymax = 1.;
  int nxmax = 8000000, nymax = 8;
  double integral = 0;
  #pragma omp parallel{
    #pragma omp single{
      for (int if = 0; if < count; if++) {
        /* code */
      }
      #pragma omp task{
        integral +=integrate(hxmin, hxmax, hymin, hymax, nxmax, nymax);
      }
    }
  }
  integral = integrate(hxmin, hxmax, hymin, hymax, nxmax, nymax);
  printf("%f\n", integral);


  tstart = clock();
  time1 += clock() - tstart;
  time1 = time1/CLOCKS_PER_SEC;
  //printf("%lf\n",time1 );

  #pragma omp parallel for
  for (size_t i = 0; i < 1000; i++) {
    //printf("I'm thread %i of %d\n", omp_get_thread_num(), omp_get_num_threads());
  }
  /* Fork a team of threads giving them their own copies of variables */
  #pragma omp parallel private(nthreads, tid)
    {
    tid = omp_get_thread_num();
    //printf("Hello World from thread = %d\n", tid);

    if (tid == 0)
      {
      nthreads = omp_get_num_threads();
      //printf("Number of threads = %d\n", nthreads);
    }
  }
}
