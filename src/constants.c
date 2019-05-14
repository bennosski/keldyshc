#include "constants.h"
#include <stdlib.h>

const int MEM_DATA_ALIGN = 16;

inline double drand()
{
  return (double)rand()/RAND_MAX;
}

inline cdouble zrand()
{
  return drand() + I*drand();
}
