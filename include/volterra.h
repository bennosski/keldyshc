#ifndef VOLTERRA_H
#define VOLTERRA_H

#include "constants.h"

void init_rcorr(cdouble * restrict rcorr);

void init_gregory_matrix_M(cdouble * restrict gmM);

void init_gregory_matrix_R(cdouble * restrict gmR);

#endif
