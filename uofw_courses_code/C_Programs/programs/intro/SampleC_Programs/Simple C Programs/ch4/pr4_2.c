/* Program 4.2 : Add three numbers */

#include <stdio.h>

void main(void)
{
  int next, sum, count;
  sum = 0;
  for (count = 1; count <= 3; count++)
  {
    scanf("%i", &next);
    sum = sum + next;
  }
  printf("%i\n", sum);
}

