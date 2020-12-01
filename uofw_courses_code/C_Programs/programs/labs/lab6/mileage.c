/* Author: Malachi Griffith
*  Date: Oct. 10 2002
*  Purpose: Calculate mileage for the user using a function mileage
*/

#include <stdio.h>

/* Function prototype */
double mileage(double miles, double gallons);

main()
{
	double miles;
	double gallons;
	double mpg;

	printf("\n\nPlease enter the miles you travelled > ");
	scanf("%lf", &miles);

	printf("Please enter the gallons of gas used > ");
	scanf("%lf", &gallons);

	mpg = mileage(miles, gallons);

	printf("The mileage of your car is: %.1f miles/gallon\n\n", mpg);	
	
	return(0);
}

double mileage(double miles, double gallons)
{
	double mpg;

	mpg = miles / gallons;
	
	return(mpg);
}
