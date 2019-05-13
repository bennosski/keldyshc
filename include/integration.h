#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "constants.h"
#include "matsubara.h"

typedef struct
{
  double * rcorr;
  double * gregory_matrix_M;
  double * gregory_matrix_R;
} integrator;

void init_integrator(integrator * restrict integ);

void MxM(const integrator * restrict integ, const matsubara * restrict A, cdouble * Cmk);

#endif


