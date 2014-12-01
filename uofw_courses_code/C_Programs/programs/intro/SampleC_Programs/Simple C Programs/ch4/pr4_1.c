/* Program 4.1 : Add three numbers */

#include <stdio.h>

void main(void)
{
  int next, sum;
  sum = 0;

  scanf("%i", &next);
  sum = sum + next;

  scanf("%i", &next);
  sum = sum + next;

  scanf("%i", &next);
  sum = sum + next;

  printf("%i\n", sum);

}

