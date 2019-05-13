#include "integration.h"
#include "params.h"
#include "matsubara.h"
A
//----------------------------------------------------
void init_integrator(integrator * restrict integ)
{
  integ->rcorr = (double *)MKL_malloc(order*order*order*sizeof(double), MEM_DATA_ALIGN);
  init_rcorr(integ->rcorr);

  integ->gregory_matrix_M = (double *)MKL_malloc(ntau*ntau*sizeof(double), MEM_DATA_ALIGN);
  init_gregory_matrix_M(integ->gregory_matrix_M);

  integ->gregory_matrix_R = (double *)MKL_malloc(nt*nt*sizeof(double), MEM_DATA_ALIGN);
  init_gregory_matrix_R(integ->gregory_matrix_R);
}

//---------------------------------------------------
void MxM(const integrator * restrict integ, const matsubara * restrict A, cdouble * Cmk)
{
  
  int iorb, korb, m, k, l;

  const int N1 = ntau;
  const int N2 = ntau*norb;
  const int N3 = ntau*norb*ntau;
  const int O1 = order;
  const int O2 = order*order;

  for(iorb=0; iorb<norb; iorb++)
    for(korb=0; korb<norb; korb++)
    {
      for(k=ntau-order; k<ntau; k++)
	for(l=0; l<order; l++)  
	  for(m=ntau-order; m<ntau-1; m++) 
	    Cmk[m + iorb*N1 + k*N2 + korb*N3] += -I*dtau * integ->rcorr[ntau-1-m + l*O1 + (ntau-1-k)*O2] * A->sig * A->M[ntau-1-l + iorb*N1 + korb*N2];

      for(m=0; m<ntau-order; m++)	 
	for(k=0; k<ntau-m; k++)
	  Cmk[m + iorb*N1 + (m+k)*N2 + korb*N3] += -I*dtau * integ->gregory_matrix_M[ntau-1-m + k*N1] * A->sig * A->M[ntau-1-k + iorb*N1 + korb*N2];
	
      for(k=0; k<order; k++)
	for(l=0; l<order; l++)
	  for(m=1; m<order; m++)
	    Cmk[m + iorb*N1 + k*N2 + korb*N3] += -I*dtau * integ->rcorr[m + l*O1 + k*O2] * A->M[l + iorb*N1 + korb*N2];

      for(m=order; m<ntau; m++)
	for(k=0; k<m+1; k++)
	  Cmk[m + iorb*N1 + k*N2 + korb*N3] += -I*dtau * integ->gregory_matrix_M[m + k*N1] * A->M[m-k + iorb*N1 + korb*N2];
    }  
}

