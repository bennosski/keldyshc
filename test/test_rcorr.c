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

void test_rcorr()
{
  printf("\n---------------------------------\n");
  printf("test rcorr weights\n---------------------------------\n");

  cdouble * rcorr = (cdouble *)MKL_malloc(order*order*order*sizeof(cdouble), MEM_DATA_ALIGN);

  char dsetname[32];
  
  sprintf(dsetname, "/rcorr_p%d", order);

  printf("dset name = %s\n", dsetname);

  dload("weights.h5", dsetname, rcorr, order*order*order);

  zsave("../test/results/test_rcorr.h5", "/rcorr", rcorr, order*order*order);

  system("exec python3 ../test/test_rcorr.py");
}




