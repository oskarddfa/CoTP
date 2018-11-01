#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>

double function(double x, double y){
  return x * x * ( - (y * y));
}

double integrate(double hxmin, double hxmax, double hymin, double hymax, int nxmax, int nymax){
  double deltahx, deltahy;
  double integral = 0;
  double nystart;
  deltahx = (hxmax - hxmin) / nxmax;
  deltahy = (hymax - hymin) / nymax;
  omp_get_num_threads();
  nystart = omp_get_thread_num() * nymax / omp_get_num_threads();
  for (int ny = nystart; ny < nymax; ny++) {
    double c = 0;
    printf("nyyyy =%d\n", ny);
    for (int nx = 0; nx < nxmax; nx++) {
      c += deltahx * deltahy / 4 * (
        function(hxmin + deltahx * nx      , hymin + deltahy * ny      ) +
        function(hxmin + deltahx * (nx + 1), hymin + deltahy * ny      ) +
        function(hxmin + deltahx * nx      , hymin + deltahy * (ny + 1)) +
        function(hxmin + deltahx * (nx + 1), hymin + deltahy * (ny + 1)));
      printf("c = %f\n", c);
    integral +=c;
    printf("nystart = %f\n", omp_get_thread_num() * nymax / omp_get_num_threads());
    //printf("int = %f\n", integral);

  return integral;

    }
    //integral +=c;
  }

  return integral;
}
int main(int argc, char *argv[]){
  double time1=0.0, tstart;
  double hxmin = 0., hymin = 0., hxmax = 1., hymax = 1.;
  int nxmax = 8000000, nymax = 8;
  double integral = 0;
  tstart = clock();
  #pragma omp parallel
  {
    #pragma omp single
    {

      #pragma omp task
      {
        integral +=integrate(hxmin, hxmax, hymin, hymax, nxmax, nymax);
        printf("thread nr:%d\n", omp_get_thread_num());
      }
    }
  }
  time1 += clock() - tstart;
  time1 = time1/CLOCKS_PER_SEC;
  printf("%f\n", time1);
  printf("%f\n", integral);


  tstart = clock();
  time1 += clock() - tstart;
  time1 = time1/CLOCKS_PER_SEC;
  //printf("%lf\n",time1 );


}
