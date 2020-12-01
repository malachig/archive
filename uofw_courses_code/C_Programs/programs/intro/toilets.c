/*toilets.c*/
/* Author: Malachi Griffith
*  Date: Sept. 20 2002 
*  Purpose: Illustrate the cost and potential savings of installing a new system of 
*	    low flush toilets in a city of a given population.
*/

#include <stdio.h>

main()
{

	int 	population, /* There are 3 people per toilet on average */
		total_toilets,
		flush_old = 15,
		flush_new = 2,
		flushes_per_day = 14,
		total_volume_old, /* volume used per day with old toilets */
		total_volume_new, /* volume used per day with new toilets */
		install_cost = 150,
		total_install_cost;


	printf("Enter the population of your city > ");
	scanf(" %d", &population);

	total_toilets = (population / 3);
	total_volume_old = total_toilets * flushes_per_day * flush_old;	
	total_volume_new = total_toilets * flushes_per_day * flush_new;
	total_install_cost = total_toilets * install_cost;

	printf("\nThe volume of water used per day with regular toilets is: %d", 
	       total_volume_old);
	printf("\nThe volume of water used per day with low usage toilets is: %d",
	       total_volume_new);
	printf("\nThe total cost to install these new toilets is: $%d\n\n",
	       total_install_cost);

}
