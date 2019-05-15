#include "volterra.h"
#include "constants.h"
#include "params.h"

void init_rcorr(double * restrict rcorr)
{
  int m, k, l;
  for(l=0; l<order; l++)
    for(k=0; k<order; k++)
      for(m=0; m<order; m++)
	rcorr[m + order*k + order*order*l] = drand();
}

void init_gregory_matrix_M(double * restrict gmM)
{
  int i, j;

  printf("ntau is %d\n", ntau);

  for(i=0; i<ntau; i++)
    for(j=0; j<ntau; j++)
      gmM[i + ntau*j] = drand();
}

void init_gregory_matrix_R(double * restrict gmR)
{
  int i, j;
  for(j=0; j<nt; j++)
    for(i=0; i<nt; i++)
      gmR[i + nt*j] = drand();
}
