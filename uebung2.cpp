#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>

Vec v;
VecCreate(MPI Comm comm,Vec *v);
VecSetSizes(Vec v, PetscInt m, PetscInt M);
VecSetFromOptions(Vec v);
VecSet(Vec x,PetscScalar value);
VecSetValues(Vec x,PetscInt n,PetscInt *indices,PetscScalar
*values,INSERT VALUES);
VecAssemblyBegin(Vec x);
VecAssemblyEnd(Vec x);
VecView(Vec x,PetscViewer v);
VecDuplicate(Vec old,Vec *new);
VecDot(Vec x,Vec y,PetscScalar *dot);

Mat A;


EPS eps; /* eigensolver context */
Mat A; /* matrix of Ax=kx */
Vec xr, xi; /* eigenvector, x */
PetscScalar kr, ki; /* eigenvalue, k */
PetscInt j, nconv;
PetscReal error;
EPSCreate( PETSC_COMM_WORLD, &eps );
EPSSetOperators( eps, A, NULL );
EPSSetProblemType( eps, EPS_NHEP );
EPSSetFromOptions( eps );
EPSSolve( eps );
EPSGetConverged( eps, &nconv );
for (j=0; j<nconv; j++) {
 EPSGetEigenpair( eps, j, &kr, &ki, xr, xi );
 EPSComputeError( eps, j, EPS_ERROR_RELATIVE, &error );
}
EPSDestroy( &eps );
