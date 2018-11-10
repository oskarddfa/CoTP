

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
