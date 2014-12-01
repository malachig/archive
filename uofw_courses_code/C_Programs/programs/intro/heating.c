/*heating.c*/
/* Author: Malachi Griffith
*  Date: Sept. 20 2002
*  Purpose: Calculates the BTUs of heat delivered to a house given the number of 
*  gallons of fuel burned and the efficiency of the heater.
*/

#include <stdio.h>

main()
{
	double efficiency,
	       conversion,
	       gallons_burned,
	       btus_delivered;

	/* Note: 42 gallons = 5,800,000 BTUs delivered at 100% efficiency */

	printf("\n\nType in the number of gallons of fuel burned > ");
	scanf(" %lf", &gallons_burned);
	printf("Type in the efficieny of your heater unit in percent> ");
	scanf(" %lf", &efficiency);
	
	conversion = 5800000.0 / 42.0; /* Number of BTUs per gallon of fuel */

	btus_delivered = (gallons_burned * conversion) * (efficiency / 100);

	printf("\nThe total BTUs delivered are %10.2f\n\n", btus_delivered);

	return(0);
}
