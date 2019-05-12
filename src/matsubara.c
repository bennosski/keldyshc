#include "matsubara.h"
#include "params.h"
#include "constants.h"

void init_matsubara(matsubara * A, int sig)
{
  A->M = (cdouble *)MKL_malloc(ntau*norb*norb*sizeof(cdouble), MEM_DATA_ALIGN);

  A->sig = sig;
}


