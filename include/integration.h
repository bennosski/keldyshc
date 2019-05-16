#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "constants.h"
#include "matsubara.h"

typedef struct
{
  cdouble * rcorr;
  cdouble * gregory_matrix_M;
  cdouble * gregory_matrix_R;
} integrator;

void init_integrator(integrator * restrict integ);

void prep_MxM(const integrator * restrict integ, const matsubara * restrict A, cdouble * Cmk);

void MxM(const matsubara * restrict A, const matsubara * restrict B, cdouble * restrict ret);

#endif


