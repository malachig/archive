/* Program 6.9 : Tolerance 4 */

#include <stdio.h>
#include <math.h>
#define STANDARD 6.37

void main(void)
{
  float length;
  int above, below, within;

  within = 0; above = 0;
  below = 0; scanf("%f", &length);

  while (length >= 0)
  {

    if (fabs(length - STANDARD) < 0.1)
      within++;
    else
      if ( (length - STANDARD) >= 0.1)
        above++;
      else
        below++;

    scanf("%f", &length);

  }

  printf("%i values are within tolerance.\n", within);
  printf("%i values are too large.\n", above);
  printf("%i values are too small.\n", below);

}

