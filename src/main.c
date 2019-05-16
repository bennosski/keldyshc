#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mkl.h>
#include <memory.h>
#include "params.h"
#include "integration.h"
#include "constants.h"
#include "util.h"
#include "matsubara.h"
#include "save.h"
#include <time.h>


//--------------------------------------------------
int main()
{
  printf("cleaning up saved files\n");
  system("exec rm results/*.h5");

  init_params();

  integrator integ;
  init_integrator(&integ);

  printf("dt %f\n", dt);
  printf("dtau %f\n", dtau);

  zsave("results/integ.h5", "/rcorr", integ.rcorr, order*order*order);
  zsave("results/integ.h5", "/gmM", integ.gregory_matrix_M, ntau*ntau);
  zsave("results/integ.h5", "/gmR", integ.gregory_matrix_R, nt*nt);

  matsubara A;
  init_matsubara(&A, -1);
  compute_G0M(&A);

  zsave("results/GM.h5", "/M", A.M, ntau*norb*norb);

  //test_dsave();
  //test_zsave();

  cdouble * Cmk = (cdouble *)MKL_calloc(ntau*norb*ntau*norb, sizeof(cdouble), MEM_DATA_ALIGN);

  clock_t time1, time2;
  time1 = clock();
  MxM(&integ, &A, Cmk);
  time2 = clock();
  printf("took : %1.4e\n", runtime(time1, time2));

  zsave("results/Cmk.h5", "/Cmk", Cmk, ntau*norb*ntau*norb);
 
  return 0;
}
