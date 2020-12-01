/* Program 5.8 : Add */

#include <stdio.h>

void main(void)
{
  float sum, next;
  sum = 0.0;
  scanf("%f", &next);
  while (next >= 0)
  {
    sum += next;
    scanf("%f", &next);
  }

  printf("Sum is %f\n", sum);
}

