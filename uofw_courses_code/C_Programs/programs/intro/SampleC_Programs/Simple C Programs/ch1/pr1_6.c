/* Program 1.6 : Area Of Circle */

#include <stdio.h>
#define PI 3.14159

void main (void)
{
  float radius;

  printf("Radius of circle?");
  scanf("%f", &radius);
  printf("The circumference is %5.3f\n", 2*PI*radius);
  printf("The area is          %5.3f\n", PI*radius*radius);

}

