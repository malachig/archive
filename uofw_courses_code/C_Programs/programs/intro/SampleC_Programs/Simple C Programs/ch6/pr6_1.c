/* Program 6.1 : Factors */

#include <stdio.h>

void main(void)
{
  int GivenInteger, i;
  scanf("%i", &GivenInteger);

  for (i=2; i<=9; i++)
    if (GivenInteger % i == 0)
      printf("%i divides into %i.\n", i, GivenInteger);

}

