/* Program 3.2 : Largest */

#include <stdio.h>

void main(void)
{
  float first, second, third, LargestSoFar, largest;
  printf("Enter three numbers:");
  scanf("%f%f%f", &first, &second, &third);

  if (first > second)
    LargestSoFar = first;
  else
    LargestSoFar = second;

  if (LargestSoFar > third)
    largest = LargestSoFar;
  else
    largest = third;

  printf("The largest number is: %.2f\n", largest);

}

