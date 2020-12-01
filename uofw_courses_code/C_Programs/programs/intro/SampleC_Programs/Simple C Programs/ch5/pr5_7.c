/* Program 5.7 : Heat Wave */

#include <stdio.h>
#define TEMP_THRESHOLD 15.0
#define SUN_THRESHOLD 10.0

void main(void)
{
  float NextTemp, NextSun;
  int day;
  FILE *WeatherFile;

  day = 0;
  WeatherFile = fopen("WEATHER.DAT", "r");
  do
  {
    day++;
    fscanf(WeatherFile, "%f%f", &NextTemp, &NextSun);
  } while (((NextTemp <= TEMP_THRESHOLD) ||
            (NextSun  <= SUN_THRESHOLD)
           ) && (day != 365) );

  if ( (NextTemp > TEMP_THRESHOLD)
    && (NextSun > SUN_THRESHOLD) )
    {
      printf("There were %i ", day-1);
      printf("days before the first good day.\n");
    }
  else
    printf("That was an exceptionally bad year.\n");

  fclose(WeatherFile);
}

