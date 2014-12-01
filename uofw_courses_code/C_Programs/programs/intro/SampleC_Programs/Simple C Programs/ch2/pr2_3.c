/* Program 2.3 : Times 2 */

#include <stdio.h>

void main(void)
{
  float dist1, dist2, dist3, mph, StartTime, TimeSoFar;

  printf("3 distances:");   scanf("%f%f%f", &dist1, &dist2, &dist3);
  printf("Start time:");    scanf("%f",&StartTime);
  printf("Average speed:"); scanf("%f", &mph);

  TimeSoFar = StartTime + dist1/mph;
  printf("Time of arrival at town 1 is %4.2f\n", TimeSoFar);

  TimeSoFar = TimeSoFar + dist2/mph;
  printf("Time of arrival at town 2 is %4.2f\n", TimeSoFar);

  TimeSoFar = TimeSoFar + dist3/mph;
  printf("Time of arrival at town 3 is %4.2f\n", TimeSoFar);

}

