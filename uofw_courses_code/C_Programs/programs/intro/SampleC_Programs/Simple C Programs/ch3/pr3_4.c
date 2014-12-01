/* Program 3.4 : Two Numbers */

#include <stdio.h>

void main(void)
{
  int larger, smaller, temporary, sum, difference;
  printf("Enter two whole numbers:");
  scanf("%i%i", &larger, &smaller);

  if (larger < smaller)
  {
    temporary = larger;
    larger = smaller;
    smaller = temporary;
  }

  sum        = larger + smaller;
  difference = larger - smaller;

  printf("\n\nThe sum is %i\n", sum);
  printf("The difference is %i\n", difference);
  printf("Numbers in order are %i, %i\n", smaller, larger);

}

