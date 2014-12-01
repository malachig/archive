/*insect.c*/
/* Author: Malachi Griffith
*  Date: Sept. 20 2002
*  Purpose: Calculate population size based on current rate of growth
*/

#include <stdio.h>

main()
{
	double population_initial,
	    population_week1,
	    population_week2;

	double growth_rate;

	printf("\n\nInput the initial population size > ");
	scanf(" %lf", &population_initial);
	printf("Input the population size after one week > ");
	scanf(" %lf", &population_week1);

	growth_rate = (double) (population_week1 / population_initial);
	printf("\n\nThe growth rate is %.2f", growth_rate);	
	population_week2 = (double) (growth_rate * population_week1);  
	printf("\nThe population after two weeks is: %.0f\n\n", population_week2);
}
