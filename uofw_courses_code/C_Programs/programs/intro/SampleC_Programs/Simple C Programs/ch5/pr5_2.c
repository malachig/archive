/* Program 5.2 : Count Days */

#include <stdio.h>
#define THRESHOLD 100.0

void main(void)
{
  float NextReading, total;
  int days;

  days = 0; total = 0;

  do
  {
    days++;
    scanf("%f", &NextReading);
    total += NextReading; 
  }
  while (total <= THRESHOLD);

  printf("Threshold total of %f has been exceeded.\n", THRESHOLD);
  printf("This occurred after %i days.\n", days);

  
}

