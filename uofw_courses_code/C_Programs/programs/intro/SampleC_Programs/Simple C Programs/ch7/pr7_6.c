/* Program 7.6 : Sort 3 numbers */

#include <stdio.h>

void order(int *a, int *b);

void main(void)
{ int first, second, third;

  scanf("%i%i%i", &first, &second, &third);
  order(&first, &second);
  order(&first, &third);
  order(&second, &third);
  printf("%6i%6i%6i\n", first, second, third);
}

void order(int *a, int *b)
{
  int temp;
  if (*a < *b)
  {
    temp = *a;
    *a = *b;
    *b = temp;
  }
}

