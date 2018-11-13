

double harmonic(double x, double dummyheight)
	{	
	return x*x;
	}


double Well(double x, double height = 20.0)
	{
	if (x<-5 || x>5){return height;}
	return 0.0;
	}
	
double Wall(double x, double height = 2.0)
	{	
	if (x>-2.0 && x<2){return height;}
	return 0.0;
	}

double Doppelmulde(double x, double a = -1.0, double b = 0.0, double c = 1.0){
  return a*x*x+b*x*x*x+c*x*x*x*x;
}

double Len_Jones(double r, double rm = 1.0, double a = 1.0){
  //Lennard Jones Potential with V(r)=4a*((rm/r)^12-2*(rm/r)^6)
  c = r * r * r * rm * rm * rm
    return 4*eps*(c*c*c*c-2*c*c*c)
    }
