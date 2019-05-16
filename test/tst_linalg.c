#include <stdio.h>
#include <complex.h>
#include "constants.h"
#define MKL_Complex16 cdouble
#include <mkl.h>
#include "save.h"


void test_linalg()
{
  printf("\nlinalg test\n");
 
  int n = 4;

  cdouble * arr1 = (cdouble *)MKL_malloc(n*n*sizeof(cdouble), MEM_DATA_ALIGN);
  cdouble * arr2 = (cdouble *)MKL_malloc(n*n*sizeof(cdouble), MEM_DATA_ALIGN);
  cdouble * ret = (cdouble *)MKL_malloc(n*n*sizeof(cdouble), MEM_DATA_ALIGN);
  
  int i,j;
  for(i=0; i<n; i++)
    for(j=0; j<n; j++)
    {
      arr1[i+n*j] = zrand();
      arr2[i+n*j] = zrand();
    }  

  cdouble alpha1=1.0, alpha2=0.0;
  cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, &alpha1, arr1, n, arr2, n, &alpha2, ret, n);

  for(i=0; i<n; i++)
  {
    for(j=0; j<n; j++)
    {
      printf("%1.2f,%1.2f ", creal(ret[i+n*j]), cimag(ret[i+n*j]));
    }
    printf("\n");
  }

  zsave("results/test_linalg.h5", "/arr1", arr1, n*n);
  zsave("results/test_linalg.h5", "/arr2", arr2, n*n);
  zsave("results/test_linalg.h5", "/ret",   ret, n*n);

  // solve A * B = C for B

  int nrhs = n;
 
  //the following solver overwrites arr1 and ret

  int * ipiv = (int *)MKL_malloc(n*sizeof(int), MEM_DATA_ALIGN);
  
  LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n, arr1, n, ipiv);

  int info;
  char trans = 'N';
  LAPACKE_zgetrs(LAPACK_COL_MAJOR, trans, n, nrhs, arr1, n, ipiv, ret, n);

  zsave("results/test_linalg.h5", "/sol", ret, n*nrhs);

  printf("test solution\n");
  
  for(i=0; i<n; i++)
  {
    for(j=0; j<n; j++)
    {
      cdouble d = ret[i+n*j]-arr2[i+n*j];
      printf("(%f,%f) ", creal(d), cimag(d));
    }
    printf("\n");
  }
	       
  MKL_free(ipiv);
  MKL_free(arr1);
  MKL_free(arr2);
  MKL_free(ret);
}
