/* Program usage: mpiexec -n 1 toy[-help] [all TAO options] */

/* ----------------------------------------------------------------------
min f=(x1-x2)^2 + (x2-2)^2 -2*x1-2*x2
s.t.     x1^2 + x2 = 2
      0 <= x1^2 - x2 <= 1
      -1 <= x1,x2 <= 2
---------------------------------------------------------------------- */

#include <petsctao.h>

static  char help[]="";

/*
   User-defined application context - contains data needed by the
   application-provided call-back routines, FormFunction(),
   FormGradient(), and FormHessian().
*/

/*
   x,d in R^n
   f in R
   bin in R^mi
   beq in R^me
   Aeq in R^(me x n)
   Ain in R^(mi x n)
   H in R^(n x n)
   min f=(1/2)*x'*H*x + d'*x
   s.t.  Aeq*x == beq
         Ain*x >= bin
*/
typedef struct {
  PetscInt n; /* Length x */
  Vec      x,xl,xu;
} AppCtx;

/* -------- User-defined Routines --------- */

PetscErrorCode InitializeProblem(AppCtx *, int xsize, int ysize);
PetscErrorCode DestroyProblem(AppCtx *);
PetscErrorCode FormFunctionGradient(Tao tao, Vec X, PetscReal *f, Vec G, void *ctx, double* bild, double alpha, int xsize);
double psf(int x, int y, double aplha);



PetscErrorCode main(int argc,char **argv)
{
  PetscErrorCode     ierr;                /* used to check for functions returning nonzeros */
  Tao                tao;
  KSP                ksp;
  PC                 pc;
  AppCtx             user;                /* application context */
  int                xsize, ysize;
  double             alpha;
  double*            bild[xsize * ysize];

  ierr = PetscInitialize(&argc,&argv,(char *)0,help);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n---- TOY Problem -----\n");CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  ierr = InitializeProblem(&user, xsize, ysize);CHKERRQ(ierr);
  ierr = TaoCreate(PETSC_COMM_WORLD,&tao);CHKERRQ(ierr);
  ierr = TaoSetType(tao,TAOIPM);CHKERRQ(ierr);
  ierr = TaoSetInitialVector(tao,user.x);CHKERRQ(ierr);
  ierr = TaoSetVariableBounds(tao,user.xl,user.xu);CHKERRQ(ierr);
  ierr = TaoSetObjectiveAndGradientRoutine(tao,FormFunctionGradient,(void*)&user);CHKERRQ(ierr);

  ierr = TaoSetFromOptions(tao);CHKERRQ(ierr);

  // ierr = TaoGetKSP(tao,&ksp);y(ierr);
  // ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  // ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
  // /*
  //     This algorithm produces matrices with zeros along the diagonal therefore we need to use
  //   SuperLU which does partial pivoting
  // */
  // ierr = PCFactorSetMatSolverType(pc,MATSOLVERSUPERLU);CHKERRQ(ierr);
  // ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  // ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /* ierr = TaoSetTolerances(tao,0,0,0);CHKERRQ(ierr); */
  ierr = TaoSolve(tao);CHKERRQ(ierr);

  ierr = DestroyProblem(&user);CHKERRQ(ierr);
  ierr = TaoDestroy(&tao);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

PetscErrorCode InitializeProblem(AppCtx *user, int xsize, int ysize)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  user->n = xsize * ysize;
  ierr = VecCreateSeq(PETSC_COMM_SELF,user->n,&user->x);CHKERRQ(ierr);
  ierr = VecDuplicate(user->x,&user->xl);CHKERRQ(ierr);
  ierr = VecDuplicate(user->x,&user->xu);CHKERRQ(ierr);
  ierr = VecSet(user->x,0.0);CHKERRQ(ierr);
  ierr = VecSet(user->xl,0.0);CHKERRQ(ierr);
  ierr = VecSet(user->xu,1.0);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode DestroyProblem(AppCtx *user)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecDestroy(&user->xl);CHKERRQ(ierr);
  ierr = VecDestroy(&user->xu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

double psf(int x,int y,double alpha){
  return exp(-1. * alpha * (x*x+y*y));
}

PetscErrorCode FormFunctionGradient(Tao tao, Vec X, PetscReal *f, Vec G, void *ctx, double* bild, double alpha, int xsize)
{
  PetscScalar       *g;
  const PetscScalar *x;
  PetscErrorCode    ierr;
  PetscFunctionBegin;

  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(G,&g);CHKERRQ(ierr);
  double f=0;
  for (size_t i = 0; i < sizeof(x) / sizeof(x[0]); i++) {
    double value=0;
    for (size_t j = 0; j < ws*ws; j++) {
      value += x[ i - ws + j%ws - ws*xsize + xsize*(j/ws)] * psf(-ws+j%ws, - ws*xsize + xsize*(j/ws), alpha);
    }
    value *= alpha / M_PI;
    f += 0.5 * (value - bild[i])*(value - bild[i]);
    g[i] = value - bild[i]
  }
  *f = f
  // *f = (x[0]-2.0)*(x[0]-2.0) + (x[1]-2.0)*(x[1]-2.0) - 2.0*(x[0]+x[1]);
  // g[0] = 2.0*(x[0]-2.0) - 2.0;
  // g[1] = 2.0*(x[1]-2.0) - 2.0;
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(G,&g);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
