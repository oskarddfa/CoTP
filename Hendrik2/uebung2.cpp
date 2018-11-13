/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-2018, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



USAGE: ./Schroedinger -n <Gridpoints> -eps_nev <#requested eigenvalues> [-Potential]



Potentiale: double harmonic(double x);
			double Well(double x, double height);
			double Wall(double x, double height);
*/


static char help[] = "Standard symmetric eigenproblem corresponding to the Laplacian operator in 1 dimension.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

#include <slepceps.h>
#include <string>
#include "Potentiale.h"
#include "ROOTplot.h"
#include <TApplication.h>
#include <iostream>




int main(int argc,char **argv)
{
  Mat            A;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki;
  Vec            xr,xi;
  PetscInt       n,i,Istart,Iend,nev,maxit,its,nconv;
  PetscErrorCode ierr;

  ierr = SlepcInitialize(&argc,&argv,(char*)0,help);if (ierr) return ierr;

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n1-D Laplacian Eigenproblem, n=%D\n\n",n);CHKERRQ(ierr);
  
  
  
  
  
  
  
  /*####################################################################
  Set parameters
  ######################################################################*/
  //x- domain
  
  double xstart = -10; double xend = 10; double height = 20.0;
  double (*V)(double,double);  // No other parameters are possible like this, i.e. Potential need to have two args input
  if (strcmp(argv[argc-1],"-Wall")==0){V = Wall;printf("POTENTIAL:\tPotentialschwelle\n");}
  else if (strcmp(argv[argc-1],"-Well")==0){V = Well;printf("POTENTIAL:\tPotentialtopf\n");}
  else if (strcmp(argv[argc-1],"-Doppelmulde")==0){V = Doppelmulde;printf("POTENTIAL:\tDoppelmulde\n");}
  else if (strcmp(argv[argc-1],"-Lennard_Jones")==0){V = Len_Jones;printf("POTENTIAL:\tLennard-Jones Potential\n");}
  else {V = harmonic;printf("POTENTIALL:\tHarmonischer Oszillator\n");}
  
  
  printf("\n");
  double h = (xend-xstart)/(n+1);
  

  /*####################################################################
  ######################################################################*/
  







  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);
  
  
  
  
  
  /*####################################################################
  Set parameters
  ######################################################################*/
  
  
  for (i=Istart;i<Iend;i++) {    
    if (i>0) { ierr = MatSetValue(A,i,i-1,-1.0/(h*h),INSERT_VALUES);CHKERRQ(ierr); }
    if (i<n-1) { ierr = MatSetValue(A,i,i+1,-1.0/(h*h),INSERT_VALUES);CHKERRQ(ierr); }
    ierr = MatSetValue(A,i,i,2.0/(h*h)+V(xstart+(i+1)*h,height),INSERT_VALUES);CHKERRQ(ierr);
  }


  //Warum in V(xstart+(i+1)*h,height) i+1 und nicht einfach i oder vielleicht i+0.5?
  
  /*####################################################################
  ######################################################################*/
  
  
  
  
  
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatCreateVecs(A,NULL,&xr);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&xi);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,A,NULL);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);
  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Get number of converged approximate eigenpairs
  */
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n");CHKERRQ(ierr);

    for (i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
      /*
         Compute the relative error associated to each eigenpair
      */
      ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);
      

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif

      if (im!=0.0) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",(double)re,(double)error);CHKERRQ(ierr);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }

  /*
     Free work space
  */
  
  
  
  /*###########################################################################
  #############################################################################
  Plotting with ROOT TGraph	*/
  


  {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      printf("\n\n\n");
      printf("Welche Eigenfunktion?\n");
      printf("k\t=\t");
      int k;
      std::cin >> k;
      if (k>nconv){ k=0;}
      // get kth converged eigenvalue eigenvector
      ierr = EPSGetEigenpair(eps,nconv-k,&kr,&ki,xr,xi);CHKERRQ(ierr);
      
      

	#if defined(PETSC_USE_COMPLEX)
      const PetscScalar *vr;
      const PetscScalar *vi;

      VecGetArrayRead(xr, &vr);
      VecGetArrayRead(xi, &vi);
	#else
      const PetscScalar *vr;
      const PetscScalar *vi;
      
      VecGetArrayRead(xr, &vr);
      VecGetArrayRead(xi, &vi);

      
	#endif
	

	auto Phisquared= new double [n];
	
	for (int i=0; i< n; i++)
		{
		Phisquared[i]= vr[i]*vr[i] + vi[i]*vi[i];
		}
	
	TApplication Dummy("App",&argc, argv);
	ROOTplot(n, Phisquared, V, xstart, xend, height);
	Dummy.Run();
  
  
	}
  
  
  
  
  
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&xr);CHKERRQ(ierr);
  ierr = VecDestroy(&xi);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
}

/*int Plothelper(int argc, char** argv) 
	{

	TApplication app("ROOT Application", &argc, argv);
	Solver(app.Argc(), app.Argv());
	app.Run();
	return 0;
}*/
