#include "matsubara.h"
#include "params.h"
#include "constants.h"

void init_matsubara(matsubara * A, int sig)
{
  A->M = (cdouble *)MKL_malloc(ntau*norb*norb*sizeof(cdouble), MEM_DATA_ALIGN);

  A->sig = sig;
}

void compute_G0M(matsubara * A)
{
  int i, a, b;
  const int N1 = ntau;
  const int N2 = ntau*norb;
  for(i=0; i<ntau; i++)
    for(a=0; a<norb; a++)
      for(b=0; b<norb; b++)
	A->M[i+N1*a+N2*b] = zrand();
}


