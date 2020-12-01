/* Author: Malachi Griffith
*  Date: Oct. 24 2002
*  Purpose: Calculate the total cost of a house for 5 years.
*/

#include <stdio.h>
#define YEARS 5

/* Functions prototypes */
double fuel_cost(double annual_fuel);
double tax(double initial_cost, double rate);
double total_cost(double initial_cost, double total_tax, 
		  double total_fuel_cost);

main()
{
	double initial;
	double tax_rate;
	double fuel;
	
	double total_taxes;
	double total_fuel;
	double grand_total;

	printf("\nEnter the initial cost of the house > ");
	scanf("%lf", &initial);
	printf("Enter the tax rate > ");
	scanf("%lf", &tax_rate);
	printf("Enter the annual fuel cost > ");
	scanf("%lf", &fuel);

	total_fuel = fuel_cost(fuel);
	total_taxes = tax(initial, tax_rate);
	grand_total = total_cost(initial, total_taxes, total_fuel);

	printf("\nThe total cost for that house is $%.2f\n\n", grand_total);
}

/*
*  Function: fuel_cost
*/
double
fuel_cost(double annual_fuel)
{
	double total;
	total = annual_fuel * YEARS;
	return(total);
}

/*
*  Function: tax
*/
double
tax(double initial_cost, double rate)
{
	double total,
	       yearly;
	yearly = initial_cost * rate;
	total = yearly * YEARS;

	return(total);
}

/*
*  Function: total_cost
*/
double total_cost(double initial_cost, double total_tax, 
		  double total_fuel_cost)
{
	double total;

	total = initial_cost + total_tax + total_fuel_cost;

	return(total);
}








