#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mkl.h>
#include <memory.h>
#include "params.h"
#include "integration.h"
#include "constants.h"
#include "matsubara.h"
#include "save.h"
#include <time.h>
#include "util.h"

void test_MxM()
{
  printf("\ntest MxM\n---------------------------------\n");

  integrator integ;

  init_integrator(&integ);

  printf("dt %f\n", dt);
  printf("dtau %f\n", dtau);

  zsave("../test/results/integ.h5", "/rcorr", integ.rcorr, order*order*order);
  zsave("../test/results/integ.h5", "/gmM", integ.gregory_matrix_M, ntau*ntau);
  zsave("../test/results/integ.h5", "/gmR", integ.gregory_matrix_R, nt*nt);

  matsubara A;
  init_matsubara(&A, -1);
  compute_G0M(&A);

  zsave("../test/results/GM.h5", "/M", A.M, ntau*norb*norb);

  cdouble * Cmk = (cdouble *)MKL_calloc(ntau*norb*ntau*norb, sizeof(cdouble), MEM_DATA_ALIGN);

  clock_t time1, time2;
  time1 = clock();
  MxM(&integ, &A, Cmk);
  time2 = clock();
  printf("c code for MxM took : %1.4e\n", runtime(time1, time2));

  zsave("../test/results/Cmk.h5", "/Cmk", Cmk, ntau*norb*ntau*norb);

  // finally run python script to compare...
  
  system("exec python3 ../test/test_MxM.py");
}

