#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <math.h>

double function(double x, double y){
  return exp(-x * x * y * y);
}

double integrate(double hxmin, double hxmax, double hymin, double hymax, int nxmax, int nymax){
  double deltahx, deltahy;
  double integral = 0;
  deltahx = (hxmax - hxmin) / nxmax;
  deltahy = (hymax - hymin) / nymax;
  for (int ny = 0; ny < nymax; ny++) {
    for (int nx = 0; nx < nxmax; nx++) {
      integral += deltahx * deltahy / 4 * (
        function(hxmin + deltahx * nx      , hymin + deltahy * ny      ) +
        function(hxmin + deltahx * (nx + 1), hymin + deltahy * ny      ) +
        function(hxmin + deltahx * nx      , hymin + deltahy * (ny + 1)) +
        function(hxmin + deltahx * (nx + 1), hymin + deltahy * (ny + 1)));
    
    }

  }

  return integral;
}


int main(int argc, char *argv[]){

  double time1=0.0, tstart;
  double hxmin = 0., hymin = 0., hxmax = 1., hymax = 1.;
  int nxmax = 8000, nymax = 8000;
  double integral = 0;
  tstart = clock();
  
  
  omp_set_dynamic(0);     // Explicitly disable dynamic teams
  omp_set_num_threads(omp_get_max_threads()); // Use n threads for all consecutive parallel regions
  //omp_set_num_threads(4); // Use n threads for all consecutive parallel regions

  int NTHREADS;
  #pragma omp parallel shared (NTHREADS, hxmin, hymin, hxmax, hymax,nxmax,nymax, integral)
  {
  #pragma omp single 
  {NTHREADS = omp_get_num_threads();printf("Currently using %d threads \n\n", NTHREADS);}
   double lochxmin, lochxmax;
      lochxmin = hxmin + omp_get_thread_num() * (hxmax - hxmin) / NTHREADS;
      lochxmax = hxmin + (omp_get_thread_num()+1) * (hxmax - hxmin) / NTHREADS;
      
     
      double partialintegral = integrate(lochxmin, lochxmax, hymin, hymax, nxmax / NTHREADS, nymax);
      printf("thread nr:%d\n", omp_get_thread_num());
      printf("\tpartialINTEGRAL: %f \n", partialintegral);
      
      #pragma omp atomic
	  integral +=partialintegral;


    printf("task ended\n" );
  }
  time1 += clock() - tstart;
  time1 = time1/CLOCKS_PER_SEC;
  
  printf("\n\n");
  printf("Time [sec]: %f\n", time1/NTHREADS);
  printf("total processorTime [sec]: %f\n", time1);
  printf("totalINTEGRAL: %f\n", integral);


  return 0;

}
