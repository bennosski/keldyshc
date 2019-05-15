#ifndef SAVE_H
#define SAVE_H

#include "constants.h"

void dsave(const char * filename, const char * dsetname, const double * restrict dset_data, int len);

void zsave(const char * filename, const char * dsetname, const cdouble * restrict dset_data, int len);

#endif
