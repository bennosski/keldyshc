

#include <stdio.h>

extern int x;

int x = 2;

void test()
{
  printf("%d\n", x);
}

void main()
{

  int a = 3;
  int x = a/4;

  printf("hello world\n");

  test();

}
