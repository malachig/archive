/* Program 3.1 : Tolerance */

#include <stdio.h>
#include <math.h>

void main(void)
{
  float standard, next; int NumberClose;
  NumberClose = 0;
  printf("Standard value?");  scanf("%f", &standard);

  printf("next test value?");  scanf("%f", &next);
  if (fabs(standard - next) < 0.1) NumberClose++;

  printf("next test value?");  scanf("%f", &next);
  if (fabs(standard - next) < 0.1) NumberClose++;

  printf("next test value?");  scanf("%f", &next);
  if (fabs(standard - next) < 0.1) NumberClose++;

  printf("%i", NumberClose);
  printf(" values are near the standard.\n");

}

