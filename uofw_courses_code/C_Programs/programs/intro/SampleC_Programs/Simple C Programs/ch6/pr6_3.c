/* Program 6.3 : Tolerance 3 */

#include <stdio.h>
#include <math.h>
#define STANDARD 6.37

void main(void)
{
  float length;
  int within, without;

  within = 0; without = 0;

  scanf("%f", &length);

  while (length >= 0)
  {
    if (fabs(length - STANDARD) < 0.1)
      within++;
    else
      without++;
    scanf("%f", &length);
  }

  printf("%i values within tolerance.\n", within);
  printf("%i values outside tolerance.\n", without);

}

