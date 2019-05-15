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

//--------------------------------------------------
void test_dsave()
{
  double arr[] = {3.0, 4.0, 5.0};
  printf("%f %f %f\n", arr[0], arr[1], arr[2]);
  char * filename = "results/test_dsave.h5";
  char * dsetname = "/dset";
  dsave(filename, dsetname, arr, 3);
}

//--------------------------------------------------
void test_zsave()
{
  cdouble arr[] = {3.0, 4.0*I, 5.0};
  printf("%f %f %f\n", arr[0], arr[1], arr[2]);
  zsave("results/test_zsave.h5", "/dset1", arr, 3);
  zsave("results/test_zsave.h5", "/dset2", arr, 3);
}

//--------------------------------------------------
double runtime(clock_t start, clock_t end)
{
  return ((double) (end-start)) / CLOCKS_PER_SEC;
}

//--------------------------------------------------
int main()
{
  printf("cleaning up saved files\n");
  /*
  int status;
  status = remove("results/test_zsave.h5");
  status = remove("results/test_dsave.h5");
  status = remove("results/GM.h5");  
  status = remove("results/integ.h5");
  status = remove("results/Cmk.h5"); 		 
  */
  system("exec rm results/*.h5");

  init_params();

  integrator integ;
  printf("success\n");

  init_integrator(&integ);

  printf("dt %f\n", dt);
  printf("dtau %f\n", dtau);

  dsave("results/integ.h5", "/rcorr", integ.rcorr, order*order*order);
  dsave("results/integ.h5", "/gmM", integ.gregory_matrix_M, ntau*ntau);
  dsave("results/integ.h5", "/gmR", integ.gregory_matrix_R, nt*nt);

  matsubara A;
  init_matsubara(&A, -1);
  compute_G0M(&A);

  zsave("results/GM.h5", "/M", A.M, ntau*norb*norb);

  test_dsave();
  test_zsave();

  cdouble * Cmk = (cdouble *)MKL_calloc(ntau*norb*ntau*norb, sizeof(cdouble), MEM_DATA_ALIGN);

  clock_t time1, time2;
  time1 = clock();
  MxM(&integ, &A, Cmk);
  time2 = clock();
  printf("took : %1.4e\n", runtime(time1, time2));

  zsave("results/Cmk.h5", "/Cmk", Cmk, ntau*norb*ntau*norb);
 
  return 0;
}
