/* Program 7.3 : Powers */

#include <stdio.h>
#include <math.h>
#define pi 3.14159
#define e 2.71828

void PrintPowers(float x, int p1, int p2);

void main(void)
{
  PrintPowers(pi, 3,6);
  PrintPowers(e,  2,7);
}

void PrintPowers(float x, int p1, int p2)
{ int p;
  printf("\n\nPowers of %.5f\n\n  n   pow(%.5f,n)\n",x,x);
  for (p=p1; p<=p2; p++)
    printf("%3i     %10.5f\n", p, pow(x,p));

}

