#ifndef VOLTERRA_H
#define VOLTERRA_H

#include "constants.h"

void init_rcorr(cdouble * restrict rcorr);

void init_gregory_matrix_M(cdouble * restrict G, int n);

#endif
