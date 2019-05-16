#include <stdio.h>
#include <mkl.h>
#include "params.h"
#include "integration.h"
#include "constants.h"
#include "matsubara.h"
#include "save.h"
#include <time.h>
#include "util.h"

void test_linalg();
void test_MxM();
void test_rcorr();

//------------------------------------------------
void test_dsave()
{
  double arr[] = {3.0, 4.0, 5.0};
  printf("%f %f %f\n", arr[0], arr[1], arr[2]);
  char * filename = "../test/results/test_dsave.h5";
  char * dsetname = "/dset";
  dsave(filename, dsetname, arr, 3);
}

//------------------------------------------------
void test_zsave()
{
  cdouble arr[] = {3.0, 4.0*I, 5.0};
  printf("%f %f %f\n", arr[0], arr[1], arr[2]);
  zsave("../test/results/test_zsave.h5", "/dset1", arr, 3);
  zsave("../test/results/test_zsave.h5", "/dset2", arr, 3);
}


void main()
{
  printf("Testing\n");

  printf("Clear results/ folder");
  system("exec rm ../test/results/*.h5");

  
  init_params();

  test_linalg();

  test_MxM();

  test_rcorr();
  
  //test_dsave();
  //test_zsave();

}
