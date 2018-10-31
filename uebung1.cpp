#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>

double function(double x, double y){
  return pow(x, 2.0) * pow(M_E, - pow(y, 2));
}

double integrate(double hxmin, double hxmax, double hymin, double hymax, int nxmax, int nymax){
  double deltahx, deltahy;
  double integral = 0;
  double time1=0.0, tstart;
  deltahx = (hxmax - hxmin) / nxmax;
  deltahy = (hymax - hymin) / nymax;
  tstart = clock();

  #pragma omp parallel for
  for (size_t nx = 0; nx < nxmax; nx++) {
    double c = 0;
    for (size_t ny = 0; ny < nymax; ny++) {
      c += deltahx * deltahy / 4 * (
        function(hxmin + deltahx * nx      , hymin + deltahy * ny      ) +
        function(hxmin + deltahx * (nx + 1), hymin + deltahy * ny      ) +
        function(hxmin + deltahx * nx      , hymin + deltahy * (ny + 1)) +
        function(hxmin + deltahx * (nx + 1), hymin + deltahy * (ny + 1)));
      //printf("c = %f\n", c);
      //printf("int = %f\n", integral);

    }
    integral +=c;
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
  int nxmax = 4000, nymax = 4000;
  double integral;
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
