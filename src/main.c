#include <stdio.h>
#include <stdlib.h>
#include <mkl.h>
#include <memory.h>
#include "params.h"
#include "integration.h"
#include "constants.h"
#include "matsubara.h"
#include "save.h"

void test(matsubara * A){
  printf("test params %d %d %f\n", norb, ntau, dtau);

  printf("%f %f %f \n", A->M[0], A->M[1], A->M[2]);
};


//--------------------------------------------------
int main()
{
  printf("hello world\n");

  printf("dt %f\n", dt);

  init_params();

  printf("dt %f\n", dt);

  integrator integ;
  init_integrator(&integ);

  matsubara A;
  init_matsubara(&A, -1);
  
  printf("norb %d\n", norb);

  A.M[0] = 0.5 + 0.3*I;
  A.M[2] = 0.2 + 0.2*I;

  test(&A);

  cdouble * Cmk = (cdouble *)MKL_calloc(ntau*norb*ntau*norb, sizeof(cdouble), MEM_DATA_ALIGN);

  MxM(&integ, &A, Cmk);

  
  return 0;
}
