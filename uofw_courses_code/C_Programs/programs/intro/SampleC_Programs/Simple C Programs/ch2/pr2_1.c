/* Program 2.1 : RoomSize3 */

#include <stdio.h>

void main (void)
{
  float length, width, height, perimeter, FloorArea;

  printf("Type in length, width and height of room:\n");
  scanf("%f%f%f", &length, &width, &height);

  perimeter = 2*(length+width);
  FloorArea = length*width;
  printf("Skirting board: %.2f m.\n", perimeter);
  printf("Wall paper: %.2f sq.m.\n", perimeter*height);
  printf("Carpet: %.2f sq.m.\n", FloorArea);
  printf("Fresh air: %.2f cu.m.\n", FloorArea*height);

}

