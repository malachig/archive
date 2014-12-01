/* Program 5.11 : Heat Wave 2 */

#include <stdio.h>
#define FALSE 0
#define TEMP_THRESHOLD 15.0
#define SUN_THRESHOLD 10.0

void main(void)
{
  float NextTemp, NextSun;
  int day;
  int WarmDayFound;

  FILE *WeatherFile;

  day = 0;
  WeatherFile = fopen("WEATHER.DAT", "r");
  WarmDayFound = FALSE;
  do {
    day++;
    fscanf(WeatherFile, "%f%f", &NextTemp, &NextSun);
    WarmDayFound =
       (NextTemp > TEMP_THRESHOLD) &&
       (NextSun  > SUN_THRESHOLD);
  } while ( (day<365) && !WarmDayFound );

  if ( WarmDayFound )
    {
      printf("There were %i ", day-1);
      printf("days before the first good day.\n");
    }
  else
    printf("That was an exceptionally bad year.\n");

  close(WeatherFile);
}

