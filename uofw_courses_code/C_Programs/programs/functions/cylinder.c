/* Author: Malachi Griffith
*  Date: Oct. 24 2002
*  Purpose: Calculate the cost of producing cylinders of different sizes
*/

#include <stdio.h>

/* Function prototype */
double surface_area(double radius, double height);
double total_cost(double surface_a, double cost_per_area, int amount);

main()
{
	double r,
	       h,
	       cost;
	int number;
	double area;
	double grand_total;

	printf("\nEnter the radius of the cylinder > ");
	scanf("%lf", &r);
	printf("Enter the height of the cylinder > ");
	scanf("%lf", &h);
	printf("Enter the cost per unit area for the material > ");
	scanf("%lf", &cost);
	printf("Enter the number of them to be produced > ");
	scanf("%d", &number);

	area = surface_area(r, h);	       
	grand_total = total_cost(area, cost, number);

	printf("The grand Total cost of the production is $%.2f\n\n"
		, grand_total);
}

/*
*  Function surface_area
*/
double
surface_area(double radius, double height)
{
	double a;
	a = 2 * 3.14159 * radius * height;
	return (a);
}

/*
*  Function total_cost
*/
double 
total_cost(double surface_a, double cost_per_area, int amount)
{
	double total;

	total = surface_a * cost_per_area * amount;
	return(total);
}


