#include <stdio.h>
#include <complex.h>
#include <mkl.h>
#include <memory.h>

typedef struct
{
  double complex * M;
  int norb;
}
matsubara;

void function(int order, double complex * arr, matsubara * A)
{
  printf("order %d\n", order);

  int i;
  for(i=0; i<4; i++)
  {
    printf("%f, %f\n", creal(A->M[i]), cimag(A->M[i]));
  }
}

void main()
{
  printf("hello world\n");

  matsubara A;
  

  A.norb = 4;

  double complex * c = (double complex *)MKL_calloc(4, sizeof(double complex *), 64);

  A.M = c;

  int x = 3;
  
  printf("hi\n");

  function(x, c, &A);

}


/*
typedef struct
{
  int * arr;
  int x; 
  double complex * arrc;
}
object;

void function(int order, object * A)
{
  printf("order %d\n", order);
  printf("A->x %d\n", A->x);
  int i;
  for(i=0; i<A->x; i++)
  {
     //printf("%d   (%f,%f)", A->arr[i], creal(A->arrc[i]), cimag(A->arrc[i]));
     printf("%d\n", A->arr[i]);
  }
}


int main()
{
  printf("hello world\n");
  int order = 6;

  object * A;

  printf("sizeof dc %d\n", sizeof(int));
  printf("sizeof dc %d\n", sizeof(double));
  printf("sizeof dc %d\n", sizeof(double complex));
  printf("sizeof dc %d\n", sizeof(*A->arrc));

  A->x = 2;
  A->arr = (int *)calloc(A->x, sizeof(int));
  //A->arrc = (double complex *)calloc(A->x, sizeof(double complex));
  
  //double complex * c = (double complex *)calloc(A->x, sizeof(double complex));

  double complex * c = (double complex *)calloc(4, sizeof(double complex));


  printf("A->x %d\n", A->x);

  function(order, A);

  return 0;
}
*/
