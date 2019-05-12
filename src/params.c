#include "params.h"

double tmax;
int nt;
double dt;
double beta;
int ntau;
double dtau;
int norb;
int order;

void init_params()
{
  tmax    = 10.0;
  nt      = 100;
  beta    = 4.0;
  ntau    = 100;
  norb    = 4;
  order   = 6;

  dt = tmax/(nt-1);
  dtau = beta/(ntau-1);
}
