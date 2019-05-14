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

void init(integrator * ptr)
{
  ptr->rcorr = (double *)MKL_malloc(order*order*order*sizeof(double), MEM_DATA_ALIGN);

  ptr->gregory_matrix_M = (double *)MKL_malloc(ntau*ntau*sizeof(double), MEM_DATA_ALIGN);

  ptr->gregory_matrix_R = (double *)MKL_malloc(nt*nt*sizeof(double), MEM_DATA_ALIGN);

  /*
  int m, k, l;
  for(l=0; l<order; l++)
    for(k=0; k<order; k++)
      for(m=0; m<order; m++)
	ptr->rcorr[m + order*k + order*order*l] = drand();
  */
  
  init_rcorr(ptr->rcorr);

  init_gregory_matrix_M(ptr->gregory_matrix_M);
}

void local_init_integrator(integrator * restrict integ)
{
  integ->rcorr = (double *)MKL_malloc(order*order*order*sizeof(double), MEM_DATA_ALIGN);

  integ->gregory_matrix_M = (double *)MKL_malloc(ntau*ntau*sizeof(double), MEM_DATA_ALIGN);

  integ->gregory_matrix_R = (double *)MKL_malloc(nt*nt*sizeof(double), MEM_DATA_ALIGN);

  init_rcorr(integ->rcorr);

  init_gregory_matrix_M(integ->gregory_matrix_M);

  init_gregory_matrix_R(integ->gregory_matrix_R);
}


//--------------------------------------------------
int main()
{
  printf("cleaning up saved files\n");
  int status;
  status = remove("results/test_zsave.h5");
  status = remove("results/test_dsave.h5");
  status = remove("results/GM.h5");  
  status = remove("results/integ.h5");
  status = remove("results/Cmk.h5"); 		 
  
  init_params();

  integrator integ;
  //init_integrator(&integ);
  //local_init_integrator(&integ);  

  printf("success\n");

  init_integrator(&integ);

  //integrator * ptr = &integ;

  //init(ptr);
  
  /*
  ptr->rcorr = (double *)MKL_malloc(order*order*order*sizeof(double), MEM_DATA_ALIGN);

  ptr->gregory_matrix_M = (double *)MKL_malloc(ntau*ntau*sizeof(double), MEM_DATA_ALIGN);

  ptr->gregory_matrix_R = (double *)MKL_malloc(nt*nt*sizeof(double), MEM_DATA_ALIGN);

  init_rcorr(ptr->rcorr);

  init_gregory_matrix_M(ptr->gregory_matrix_M);
  */

  printf("dt %f\n", dt);
  printf("dtau %f\n", dtau);

  /*
  integrator integ;
  init_integrator(&integ);
  */

  printf("h1\n");

  return 0;

  dsave("results/integ.h5", "/rcorr", integ.rcorr, order*order*order);
  dsave("results/integ.h5", "/gmM", integ.gregory_matrix_M, ntau*ntau);
  dsave("results/integ.h5", "/gmR", integ.gregory_matrix_R, nt*nt);

  printf("h2\n");

  matsubara A;
  init_matsubara(&A, -1);
  compute_G0M(&A);

  printf("h3\n");

  zsave("results/GM.h5", "/M", A.M, ntau*norb*norb);

  printf("h4\n");

  test_dsave();
  test_zsave();

  cdouble * Cmk = (cdouble *)MKL_calloc(ntau*norb*ntau*norb, sizeof(cdouble), MEM_DATA_ALIGN);
  MxM(&integ, &A, Cmk);

  zsave("results/Cmk.h5", "/Cmk", Cmk, ntau*norb*ntau*norb);
 
  return 0;
}
