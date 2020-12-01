/* Program 2.7 : Powers */

#include <stdio.h>
#include <math.h>

void main (void)
{
  float edge, surface, volume;

  printf("Length of edge?\n");
  scanf("%f", &edge);

  surface = 6 * pow(edge, 2);
  volume  = pow(edge, 3);
  printf("Surface area is %.2f, volume is %.2f\n", surface, volume);

}

