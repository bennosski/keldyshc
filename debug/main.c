#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <mkl.h>
#include <memory.h>

int MEM_DATA_ALIGN=64;
typedef double complex cdouble;
//--------------------------------------------------
typedef struct
{
  cdouble * M;
  int sig;
  int norb;
  int ntau;
  double dtau; 
}
matsubara;

typedef struct
{
  cdouble * rcorr;
  cdouble * gregory_matrix_M;
  cdouble * gregory_matrix_R;
}

//--------------------------------------------------
void MxM(const int order, integrator * integ, matsubara * restrict A, cdouble * Cmk)
{
  double dtau = A->dtau;
  int ntau = A->ntau;
  int norb = A->norb;
  int sig = A->sig;

  int iorb, korb, m, k, l;
  for(iorb=0; iorb<norb; iorb++)
    for(korb=0; korb<norb; korb++)
    {
      for(m=ntau-order; m<ntau-1; m++)
        for(k=ntau-order; k<ntau; k++)
	  for(l=0; l<order; l++)  
	    Cmk[m iorb k korb] += I*dtau * integ->rcorr[ntau-1-m l ntau-1-k] * sig * A->M[ntau-1-l iorb korb];

	 
      for(m=0; m<ntau-order; m++)
        for(k=0; k<ntau-m; k++)
	  Cmk[m iorb m+k korb] += -I*dtau * integ->gregory_matrix_M[ntau-1-m k] * sig * A->M[ntau-1-k iorb korb];
	
      for(m=1; m<order; m++)
	for(k=0; k<order; k++)
	  for(l=0; l<order; l++)
	    Cmk[m iorb k korb] += -I*dtau * integ->rcorr[m l k] * A->M[l irob, korb];

      for(m=order; m<ntau; m++)
	for(k=0; k<m+1; k++)
	  Cmk[m iorb k korb] += -I*dtau * integ->gregory_matrix_M[m k] * A->M[m-k iorb korb];
    }

  
  /*
  int x;  
  for(x=0; x<A->norb; x++)
  {
    printf("%f %f\n", creal(A->M[x]), cimag(A->M[x]));
  }
  */
}
//--------------------------------------------------
int main()
{
  printf("hello world\n");

  matsubara A;

  A.ntau = 1000;
  A.dtau = 0.1;
  A.norb = 4;
  A.sig = -1;

  int order = 6;

  A.M = (cdouble *)MKL_malloc(A.ntau*A.norb*A.norb*sizeof(cdouble), MEM_DATA_ALIGN);

  printf("A norb %d\n", A.norb);

  A.M[0] = 0.5 + 0.3*I;
  A.M[2] = 0.2 + 0.2*I;

  cdouble * Cmk = (cdouble *)MKL_calloc(A.ntau*A.norb*A.ntau*A.norb, sizeof(cdouble), MEM_DATA_ALIGN);

  MxM(order, &A, Cmk);

  return 0;
}
