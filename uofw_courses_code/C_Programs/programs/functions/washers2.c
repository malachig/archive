/* Author: Malachi Griffith
*  Date: Sept. 8 2002
*  Purpose: Computes the weight of a batch of flat washers
*	    Also find the material needed to produce these washers, 
*	    if they are stamped out of rectangular pieces, and the material 
*	    left over. 
*/

#include <stdio.h>
#define PI 3.14159

int
main(void)
{
	double hole_diameter;	/* input diameter of hole 	*/
	double edge_diameter;	/* input - diameter of outer edge */
	double thickness;	/* input - thickness of washer */
	double density;		/* input - density of material used */
	double quantity;	/* input - number of washers made */
	double weight;		/* output - weight of washer batch */
	double hole_radius;	/* radius of hole */
	double edge_radius;	/* radius of outer edge */
	double rim_area;	/* area of rim	*/
	double unit_weight;	/* weight of 1 washer */

	double unit_square_area,
	       total_square_area,
	       square_weight,
	       wasted_weight;
 
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
	
	/* Compute the rim area. */
	hole_radius = hole_diameter / 2.0;
	edge_radius = edge_diameter / 2.0;
	rim_area = (PI * edge_radius * edge_radius) - (PI * hole_radius * hole_radius);
	
	/* Compute the weight of a flat washer */
	unit_weight = rim_area * thickness * density;
	
	/* Compute the weight of a batch of washers. */
	weight = unit_weight * quantity;

	/* Compute the amount of material needed to produce the washers. */
	unit_square_area = edge_diameter * edge_diameter;
	total_square_area = unit_square_area * quantity;
	square_weight = total_square_area * thickness * density;
	wasted_weight = square_weight - weight;

	/* Display the weight of the batch of washers. */
	printf("\nThe expected weight of the batch is %.2f", weight);
	printf(" grams.\n");

	printf("\nThe weight of squares required to produce these is : %.2f", square_weight);
	printf("\nThe weight wasted by cutting out the washers is : %.2f\n\n", wasted_weight);
	
	return(0);
}
