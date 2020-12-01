/* Program 4.5 : Commission Table */

#include <stdio.h>
#define LOWER_RATE 0.125
#define HIGHER_RATE 0.155
#define TOP_OF_LOWER_RANGE 100

void main(void)
{
  int sale;
  printf("Sale   Commission\n\n");

  for (sale = 1; sale <= TOP_OF_LOWER_RANGE; sale++)
    printf("%4i%13.2f\n", sale, sale * LOWER_RATE);

  for (sale = TOP_OF_LOWER_RANGE+1; sale <= 200; sale++)
    printf("%4i%13.2f\n", sale, sale * HIGHER_RATE);

}

