#include "volterra.h"
#include "constants.h"
#include "params.h"
#include "util.h"
#include "save.h"
#include <stdlib.h>

/*
void start_poly_w(double * wpoly)
{
  //w = np.zeros((p,p))

  if(order==2)
  {
    double w[] =
      {0.0, 0.5,
       0.0, 0.5};
    memcpy(wpoly, &w, order*order*sizeof(double));
  }
  else if(order==3)
  {
    double w[] =
      {0.0, 5.0/12.0, 1.0/3.0
       0.0, 2.0/3.0,  4.0/3.0
       0.0,-1.0/12.0, 1.0/3.0};
    memcpy(wpoly, &w, order*order*sizeof(double));
  }
  else if(order==4)
  {
    double w[] = 
      {0.0, 3.0/8.0, 1.0/3.0, 0.0, 
       1.0, 19.0/24.0, 4.0/3.0, 0.0, 
       1.0, -5.0/24.0, 1.0/3.0, 0.0, 
       1.0, 1.0/24.0, , 0.0};
    memcpy(wpoly, &w, order*order*sizeof(double));    
  }
  else if(order==5)
  {
    double w[] = 
      {0.0, 0.0, 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0, 0.0, 0.0};
    memcpy(wpoly, &w, order*order*sizeof(double));
  }
  else if(order==6)
  {
    double w[] = 
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
       1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    memcpy(wpoly, &w, order*order*sizeof(double));
  }

  
}  

  /*
  else if (p==3)
    w[1,:] = [5.0/12.0,2.0/3.0,-1.0/12.0]
    w[2,:] = [1.0/3.0,4.0/3.0,1.0/3.0]
    elif p==4:
        w[1,:] = [3.0/8.0,19.0/24.0,-5.0/24.0,1.0/24.0]
        w[2,:] = [1.0/3.0,4.0/3.0,1.0/3.0,0.0]
        w[3,:] = [3.0/8.0,9.0/8.0,9.0/8.0,3.0/8.0]
    elif p==5:
        w[1,:] = [251.0/720.0,323.0/360.0,-11.0/30.0,53.0/360.0,-19.0/720.0]
        w[2,:] = [29.0/90.0,62.0/45.0,4.0/15.0,2.0/45.0,-1.0/90.0]
        w[3,:] = [27.0/80.0, 51.0/40.0, 9.0/10.0, 21.0/40.0, -(3.0/80.0)]
        w[4,:] = [14.0/45.0, 64.0/45.0, 8.0/15.0, 64.0/45.0, 14.0/45.0]
    elif p==6:
        w[1,:] = [95.0/288.0, 1427.0/1440.0, -(133.0/240.0), 241.0/720.0, -(173.0/1440.0), 3.0/160.0]
        w[2,:] = [14.0/45.0, 43.0/30.0, 7.0/45.0, 7.0/45.0, -(1.0/15.0), 1.0/90.0]
        w[3,:] = [51.0/160.0, 219.0/160.0, 57.0/80.0, 57.0/80.0, -(21.0/160.0), 3.0/160.0]
        w[4,:] = [14.0/45.0, 64.0/45.0, 8.0/15.0, 64.0/45.0, 14.0/45.0, 0.0]
        w[5,:] = [95.0/288.0, 125.0/96.0, 125.0/144.0, 125.0/144.0, 125.0/96.0, 95.0/288.0]
    return w
  */

/*
void AdamsMoultonW(p):
    if p==1:
        w = [1.0]
        return w
    elif p==2:
        w = [0.5,0.5]
        return w
    elif p==3:
        w = [-(1.0/12.0), 2.0/3.0, 5.0/12.0]
        return w
    elif p==4:
        w = [1.0/24.0, -(5.0/24.0), 19.0/24.0, 3.0/8.0]
        return w
    elif p==5:
        w = [-(19.0/720.0), 53.0/360.0, -(11.0/30.0), 323.0/360.0, 251.0/720.0]
        return w
    elif p==6:
        w = [3.0/160.0, -(173.0/1440.0), 241.0/720.0, -(133.0/240.0), 1427.0/1440.0, 95.0/288.0]
        return w

void OmegaW(p):
    if p==1:
        w = [1.0]
    elif p==2:
        w = [0.5,1.0]
    elif p==3:
        w = [5.0/12.0, 13.0/12.0, 1.0]
    elif p==4:
        w = [3.0/8.0, 7.0/6.0, 23.0/24.0, 1.0]
    elif p==5:
        w = [251.0/720.0, 299.0/240.0, 211.0/240.0, 739.0/720.0, 1.0]
    elif p==6:
        w = [95.0/288.0, 317.0/240.0, 23.0/30.0, 793.0/720.0, 157.0/160.0, 1.0]
    return w

void weights(p):
    wpoly = StartPolyW(p)

    wam = AdamsMoultonW(p)

    omega = OmegaW(p)

    wstart = np.zeros((2*p-1,p))

    wstart[0:p,0:p] = wpoly
    for k in range(0,p-1):
        for j in range(0,k+1):
            wstart[p+k,j] = wstart[p+k-1,j]
        for j in range(k+1,p):
            wstart[p+k,j] = wstart[p+k-1,j] + wam[j-k-1]

    return wstart,omega

*/


void init_rcorr(cdouble * restrict rcorr)
{
  char dsetname[32];
  sprintf(dsetname, "/rcorr_p%d", order);
  dload("weights.h5", dsetname, rcorr, order*order*order);
}

void init_gregory_matrix_M(cdouble * restrict gmM)
{

  cdouble wstart[(2*order-1)*order];
  cdouble omega[order];

  char dsetname[32];
  sprintf(dsetname, "/wstart_p%d", order);
  dload("weights.h5", dsetname, wstart, (2*order-1)*order);
  sprintf(dsetname, "/omega_p%d", order);
  dload("weights.h5", dsetname, omega, order);

  int i, j;
  for(i=0; i<ntau; i++)
    for(j=0; j<ntau; j++)
      gmM[i + ntau*j] = drand();
}

void init_gregory_matrix_R(cdouble * restrict gmR)
{
  cdouble wstart[(2*order-1)*order];
  cdouble omega[order];

  char dsetname[32];
  sprintf(dsetname, "/wstart_p%d", order);
  dload("weights.h5", dsetname, wstart, (2*order-1)*order);
  sprintf(dsetname, "/omega_p%d", order);
  dload("weights.h5", dsetname, omega, order);

  int i, j;
  for(j=0; j<nt; j++)
    for(i=0; i<nt; i++)
      gmR[i + nt*j] = drand();
}

