/* Program 6.2 : Tolerance 2 */

#include <stdio.h>
#include <math.h>

void main(void)
{
  float standard, next;
  int NumberClose, count;
  NumberClose = 0;
  scanf("%f", &standard);

  for (count = 1; count <= 3; count++)
  {
    scanf("%f", &next);
    if (fabs(standard - next) < 0.1)
      NumberClose++;
  }

  printf("%i", NumberClose);
  printf(" values are near the standard.\n");

}

