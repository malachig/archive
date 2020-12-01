/* Program 2.6 : Roof */

#include <stdio.h>
#include <math.h>

void main (void)
{
  float pitch, length, width, AreaOfRoof;

  printf("Type in pitch of roof (degrees), length and width:\n");
  scanf("%f%f%f", &pitch, &length, &width);

  pitch = pitch * 3.14159/180;
  AreaOfRoof = width / cos(pitch) * length;

  printf("Area of roof is %f sq.m.\n", AreaOfRoof);

}

