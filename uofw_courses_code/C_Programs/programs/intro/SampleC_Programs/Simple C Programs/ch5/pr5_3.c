/* Program 5.3 : Next Train */

#include <stdio.h>

void main(void)
{
  int earliest, depart, arrive;

  printf("Earliest possible depart time?");  scanf("%i", &earliest);

  printf("Now input the timetable:\n");

  do
    scanf("%i%i", &depart, &arrive);
  while (depart < earliest);

  printf("Your train leaves at %i ", depart);
  printf("and arrives at %i\n", arrive);
}

