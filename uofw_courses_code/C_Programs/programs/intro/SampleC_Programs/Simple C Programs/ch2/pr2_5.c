/* Program 2.5 : Journey */

#include <stdio.h>

void main (void)
{
  int StartTime, hours, mins;
  float dist, mph, FloatHours;

  printf("Distance:");   scanf("%f", &dist);
  printf("Start time:");    scanf("%i",&StartTime);
  printf("Average speed:"); scanf("%f", &mph);

  hours = StartTime / 100;
  mins  = StartTime % 100;
  FloatHours = hours + mins/60.0;

     /* next statement calculates time of arrival */
     /* as a float number of hours after midnight  */
  FloatHours = FloatHours + dist/mph;

     /* time now converted back into minutes */
     /* rounded to nearest minute */
  mins  = FloatHours*60 + 0.5;

     /* now convert minutes to hours and minutes  */
  hours = mins/60;
  mins = mins % 60;
  printf("Arrival time: %2i:%02i\n", hours, mins);

}

