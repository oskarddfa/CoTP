#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>

double function(double x, double y){
  return 1;//x * x * (  (y * y));
}

double integrate(double hxmin, double hxmax, double hymin, double hymax, int nxmax, int nymax){
  double deltahx, deltahy;
  double integral = 0;
  deltahx = (hxmax - hxmin) / nxmax;
  deltahy = (hymax - hymin) / nymax;
  printf("%d, %d\n",nymax,nxmax );
  for (int ny = 0; ny < nymax; ny++) {
    for (int nx = 0; nx < nxmax; nx++) {
      integral += deltahx * deltahy / 4 * (
        function(hxmin + deltahx * nx      , hymin + deltahy * ny      ) +
        function(hxmin + deltahx * (nx + 1), hymin + deltahy * ny      ) +
        function(hxmin + deltahx * nx      , hymin + deltahy * (ny + 1)) +
        function(hxmin + deltahx * (nx + 1), hymin + deltahy * (ny + 1)));
    //printf("nystart = %f\n", omp_get_thread_num() * nymax / omp_get_num_threads());
    //printf("int = %f\n", integral);
    }
    //integral +=c;
  }

  return integral;
}
int main(int argc, char *argv[]){
  int k=0;

  double time1=0.0, tstart;
  double hxmin = 0., hymin = 0., hxmax = 1., hymax = 1.;
  int nxmax = 8000, nymax = 8000;
  double integral = 0;
  tstart = clock();
  //#pragma omp parallel
  { double lochxmin, lochxmax, lochymin, lochymax;
    int locnxmax = 32000, locnymax = 2000;
      lochxmin = hxmin;//hxmin + omp_get_thread_num() * (hxmax - hxmin) / omp_get_num_threads();
      lochxmax = hxmax;//hxmin + (omp_get_thread_num()+1) * (hxmax - hxmin) / omp_get_num_threads();
      //lochymin = hymin + omp_get_thread_num() * (hymax - hymin) / omp_get_num_threads();
      //lochymax = hymin + (omp_get_thread_num()+1) * (hymax - hymin) / omp_get_num_threads();
      lochymin = hymin;
      lochymax = hymax;
      integral += integrate(lochxmin, lochxmax, lochymin, lochymax, locnxmax, locnymax);
      printf("thread nr:%d\n", omp_get_thread_num());
      //printf("INTEGRAL:%f, %f\n", integral2, lochymax);

    #pragma omp taskwait
    printf("task ended\n" );
  }
  time1 += clock() - tstart;
  time1 = time1/CLOCKS_PER_SEC;
  printf("%f\n", time1);
  printf("%f\n", integral);


  tstart = clock();
  time1 += clock() - tstart;
  time1 = time1/CLOCKS_PER_SEC;
  //printf("%lf\n",time1 );
  scanf("%d\n", (int) k);
  return 0;

}
