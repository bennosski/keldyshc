#ifndef MATSUBARA_H
#define MATSUBARA_H

#include "constants.h"

typedef struct
{
  cdouble * M;
  int sig;
} matsubara;

void init_matsubara(matsubara * mat, int sig);

#endif


