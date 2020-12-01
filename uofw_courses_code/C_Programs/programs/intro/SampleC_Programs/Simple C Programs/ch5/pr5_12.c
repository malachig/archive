/* Program 5.12 : Next Train 3*/

#include <stdio.h>
#define FALSE 0


void main(void)
{
  int earliest, depart, arrive;
  int TrainFound;
  FILE *timetable;

  printf("Earliest possible depart time?");  scanf("%i", &earliest);

  timetable = fopen("TIMES.DAT", "r");

  TrainFound=FALSE;
  do {
    fscanf(timetable, "%i%i", &depart, &arrive);
    if (!feof(timetable)) TrainFound = (depart >= earliest);

  } while (!feof(timetable) && !TrainFound);

  if (TrainFound)
  {
    printf("Your train leaves at %i ", depart);
    printf("and arrives at %i\n", arrive);
  }
  else
    printf("No train available.\n");

  fclose(timetable);

  
}

