/* population.c */

/* Author: Malachi Griffith
*  Date: Oct. 5 2002
*  Purpose: Use a loop to determine when a city's population exceeds
*           30000, at a growth rate of 10% per year.
*	ALSO INCLUDES A SIMPLE FUNCTION FOR ROUNDING NUMBERS
*/

#include <stdio.h>
#define round(x) ((x)>=0?(long) ((x)+0.5):(long)((x)-0.5))
main()
{
	double init_pop = 9870.0;
	double max_pop = 30000.0;
	double current;	
	double growth_pop = 0.0;	
	int count_years = 0;
	double growth_rate = 0.10;
	double y;

	for (current = init_pop;
	     current <= max_pop;
	     current += growth_pop)
	{
	growth_pop = round((current * 0.10));		
	printf("\nPopulation = %.0f after year %d", current, count_years);
	count_years++;
	}
	printf("\nThe population exceeded 30000 people after %d years\n",
	       count_years);
	

	y = round(1.1);
	printf ("\n\nRounding of 1.1 results in: %.1f", y);
	y = round(1.49);
	printf ("\nRounding of 1.49 results in: %.1f", y);
	y = round(1.51);
	printf ("\nRounding of 1.51 results in: %.1f", y);
	y = round(1.71);
	printf ("\nRounding of 1.71 results in: %.1f\n\n", y);


}



