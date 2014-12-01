/* Program 5.4 : Next Train 2 */

#include <stdio.h>

void main(void)
{
  int earliest, depart, arrive;
  FILE *timetable;

  printf("Earliest possible depart time?");  scanf("%i", &earliest);

  timetable = fopen("TIMES.DAT", "r");

  do
    fscanf(timetable, "%i%i", &depart, &arrive);
  while (depart < earliest);

  printf("Your train leaves at %i ", depart);
  printf("and arrives at %i\n", arrive);

  fclose(timetable);

  
}

