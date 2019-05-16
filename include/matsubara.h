#ifndef MATSUBARA_H
#define MATSUBARA_H

#include "constants.h"

typedef struct
{
  cdouble * M;
  cdouble * delta;
  int sig;
} matsubara;

void init_matsubara(const matsubara * mat, int sig);

#endif


