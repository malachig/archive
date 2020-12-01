/* Program 3.5 : Two Numbers */

#include <stdio.h>

void main(void)
{
  int first, second, larger, smaller, sum, difference;
  printf("Enter two whole numbers:");
  scanf("%i%i", &first, &second);

  if (first < second)
  {
     smaller = first;
     larger  = second;
  }
  else
  {
    smaller = second;
    larger  = first;
  }

  sum        = larger + smaller;
  difference = larger - smaller;

  printf("\n\nThe sum is %i\n", sum);
  printf("The difference is %i\n", difference);
  printf("Numbers in order are %i, %i\n", smaller, larger);
}

