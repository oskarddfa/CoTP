
#include <cmath>

//Funktion die im Paper von Boys benutzt wird
double F(double x){
  return 1 / sqrt(x) * (erf(sqrt(x))- erf(0));
}

double calculate_AB_2(double Ax, double Ay, double Az, double Bx, double By, double Bz){
  double AB_2 = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);
  return AB_2;
}

// Freies Teilchen
double free(double a, double b, double AB_2){
  // Für A, B braucht man im Normalfall nur Ax und Bx zu definieren
  double result;
  // result wurde in 2 Zeilen geteilt zur Übersicht
  result = (3*a*b) / (a+b) - (2 * AB_2 * a*a * b*b) / ((a+b) * (a+b));
  result = result * pow(M_PI / (a+b), 1.5) * exp( - AB_2*a*b / (a+b));
  return result;
}

// Elektron - Kern Wechselwirkung
double nucleus(double a, double b, double AB_2, double KP_2){

  double vorfaktor, result;

  vorfaktor = ((2 * M_PI) / (a + b)) * exp(- (AB_2 * a * b) / (a+b));
  result = vorfaktor * F(KP_2 * (a + b));
  return result;
}

/*
// Elektron - Elektron Wechselwirkung
double ee_interaction(double a, double b, double c, double d, double AB_2, double CD_2){

  double result;
  result = (2 * pow(M_PI, 2.5) / ((a+b) * (c+d) * sqrt(a+b+c+d)));
  result = result * exp( - (AB_2 * a*b) / (a+b) - (CD_2 *c*d) / (c+d));
  result = result * F((PQ_2 * (a+b)*(c+d)) / (a+b+c+d));
  return result;
  }
*/


double hamiltonian(double a[], double b[], double Ax[], double Ay[], double Az[],
  double Bx[], double By[], double Bz[], double Kx[], double Ky[], double Kz[], int N_nuclei, int N_particles=1){
  /*
  a, b und A, B sind die Variablen der Wellenfunktionen
  K ist die Position der Kerne
  N_nuclei ist die Anzahl der Kerne
  N_particles ist die Anzahl der Teilchen
  Bei einem einzelnen Teilchen muessen &a, &b, &Ax usw. weitergegeben werden

  */
  double result=0.0;

  for (int i = 0; i < N_particles; i++) {
    double sum_of_nuclei = 0.0;
    // AB_2 ist sum_x,y,z  (A_i - B_i)^2
    double AB_2 = calculate_AB_2(Ax[i], Ay[i], Az[i], Bx[i], By[i], Bz[i]);
    // Px = (a * Ax + b * Bx) / (a + b) ist der Eintrag zum Ausrechnen von KP_2
    double Px = (a[i]*Ax[i]+b[i]*Bx[i]) / (a[i]+b[i]);
    double Py = (a[i]*Ay[i]+b[i]*By[i]) / (a[i]+b[i]);
    double Pz = (a[i]*Az[i]+b[i]*Bz[i]) / (a[i]+b[i]);

    for (int j = 0; j < N_nuclei; j++) {
      double KP_2 = calculate_AB_2(Kx[j], Ky[j], Kz[j], Px, Py, Pz);
      sum_of_nuclei += nucleus(a[i], b[i], AB_2, KP_2);
    }
    /*
    double sum_of_ee_interactions = 0.0;
    FUER MEHRERE TEILCHEN
    for (size_t j = 0; j < N_particles; j++) {
      double CD_2 = calculate_AB_2(Ax[j], Ay[j], Az[j], Bx[j], By[j], Bz[j]);
      ee_interaction(a[i], b[i], a[j], b[j])
    }
    */

    result += free(a[i], b[i], AB_2) + sum_of_nuclei;
  }
  return result;
}
