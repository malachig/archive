/* Program 2.2 : Times */

#include <stdio.h>

void main (void)
{
  float dist1, dist2, dist3, mph;
  float StartTime, TimeAt1, TimeAt2, TimeAt3;

  printf("3 distances:");   scanf("%f%f%f", &dist1, &dist2, &dist3);
  printf("Start time:");    scanf("%f",&StartTime);
  printf("Average speed:"); scanf("%f", &mph);

  TimeAt1 = StartTime + dist1/mph;
  printf("Time of arrival at town 1 is %4.2f\n", TimeAt1);

  TimeAt2 = StartTime + (dist1+dist2)/mph;
  printf("Time of arrival at town 2 is %4.2f\n", TimeAt2);

  TimeAt3 = StartTime + (dist1+dist2+dist3)/mph;
  printf("Time of arrival at town 3 is %4.2f\n", TimeAt3);

}

