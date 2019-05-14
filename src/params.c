#include "params.h"
#include <stdio.h>

double tmax;
int nt;
double beta;
int ntau;
int norb;
int order;

double dt;
double dtau;

void parse_line(char * line, FILE * fp)
{
  fgets(line, 256, fp);  
  int i;
  for(i=0; i<strlen(line)-1; i++)
  {
    if(line[i]==' ')
    {
      line[i] = '\0';
      return;
    }
  }  
  return;
}

void read_params()
{
  FILE *fp;
  char line[256];
  
  fp = fopen("params.txt", "r");

  parse_line(line, fp); tmax  = atof(line);
  parse_line(line, fp); nt    = atoi(line);
  parse_line(line, fp); beta  = atof(line);
  parse_line(line, fp); ntau  = atoi(line);
  parse_line(line, fp); norb  = atoi(line);
  parse_line(line, fp); order = atoi(line);

  printf("params\n");
  printf("%f %d %f %d %d %d\n", tmax, nt, beta, ntau, norb, order);

  fclose(fp);
}

void init_params()
{
  read_params();

  /*
  tmax    = 10.0;
  nt      = 20;
  beta    = 4.0;
  ntau    = 128;
  norb    = 3;
  order   = 6;
  */

  dt = tmax/(nt-1);
  dtau = beta/(ntau-1);
}
