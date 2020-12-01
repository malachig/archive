/* emissions.c */

/* Author: Malachi Griffith
*  Date: Oct. 25 2002
*  Purpose: Determine whether user has exceed allowable emissions
*/

#include <stdio.h>

main()
{
	int pollutant;
	double odometer_reading;
	double grams_per_km;

	double first_50000[4] = {3.4, 0.31, 0.4, 0.25};
	double second_50000[4] = {4.2, 0.39, 0.5, 0.31};

	printf("\n(1) Carbon monoxide");
	printf("\n(2) Hydrocarbons");
	printf("\n(3) Nitrogen oxides");
	printf("\n(4) Nonmethane hydrocarbons");

	printf("\n\nEnter the pollutant number > ");
	scanf("%d", &pollutant);
	printf("Enter number of grams emitted per km > ");
	scanf("%lf", &grams_per_km);
	printf("Enter the odometer reading > ");
	scanf("%lf", &odometer_reading);

	if (odometer_reading <= 50000)
	 	if (grams_per_km > first_50000[pollutant -1])
			{	
			printf("\nEmissions exceed permitted level");
			printf(" of %.2f grams/km\n\n",
				 first_50000[pollutant - 1]);
 			}
		else
			printf("\nYour ok dude!\n\n");	
	else
		if(grams_per_km > second_50000[pollutant - 1])
			{
			printf("\nEmissions exceed permitted level");
			printf(" of %.2f grams/km\n\n",
				second_50000[pollutant - 1]);
			}
		else
			printf("\nYour ok dude!\n\n");
}



