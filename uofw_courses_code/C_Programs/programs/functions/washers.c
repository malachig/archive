/* Author: Malachi Griffith
*  Date: Sept. 8 2002
*  Purpose: Computes the weight of a batch of flat washers
*/

#include <stdio.h>
#define PI 3.14159

/* Function Prototypes */
double find_area(double hole_diameter, double edge_diameter);
double find_unit_weight(double thickness, double rim_area, double density); 
double find_total_weight(double unit_weight, double quantity);

int
main(void)
{
	double hole_diameter;	/* input diameter of hole 	*/
	double edge_diameter;	/* input - diameter of outer edge */
	double thickness;	/* input - thickness of washer */
	double density;		/* input - density of material used */
	double quantity;	/* input - number of washers made */
	double rim_area;	/* area of rim	*/
	double weight;		/* output - weight of washer batch */
	double unit_weight;	/* weight of 1 washer */

	/* Get the inner diameter, outer diameter, and thickness. */
	printf("\n\nInner diameter in centimeters > ");
	scanf("%lf", &hole_diameter);
	printf("Outer diameter in centimeters > ");
	scanf("%lf", &edge_diameter);
	printf("Thickness in centimeters > ");
	scanf("%lf", &thickness);

	/* Get the material density and quantity manufactures. */
	printf("Material density in grams per cubis centimeter > ");
	scanf("%lf", &density);
	printf("Quantity in batch > ");
	scanf("%lf", &quantity);

	rim_area = find_area(hole_diameter, edge_diameter);	

	unit_weight = find_unit_weight(thickness, rim_area, density);

	weight = find_total_weight(unit_weight, quantity);	

	/* Display the weight of the batch of washers. */
	printf("\nThe expected weight of the batch is %.2f", weight);
	printf(" grams.\n");
	
	return(0);
}


/* Compute the rim area. */
double 
find_area(double hole_diameter, double edge_diameter)
{
	double hole_radius;	/* radius of hole */
	double edge_radius;	/* radius of outer edge */
	double rim_area;	/* area of rim	*/
	
	hole_radius = hole_diameter / 2.0;
	edge_radius = edge_diameter / 2.0;
	rim_area = (PI * edge_radius * edge_radius) - 
	           (PI * hole_radius * hole_radius);
	return(rim_area);
}


/* Compute the weight of a flat washer */
double 
find_unit_weight(double thickness, double rim_area, double density) 
{
	double unit_weight;
	unit_weight = rim_area * thickness * density;
	
	return(unit_weight);
}


/* Compute the weight of a batch of washers. */
double 
find_total_weight(double unit_weight, double quantity)
{
	double weight;	
	weight = unit_weight * quantity;

	return(weight);
}
