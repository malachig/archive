/* Author: Malachi Griffith
*  Date: Oct. 12 2002
*  Purpose: Calculates the distance between two points on a 2-d plane
*	    by using their cartesian coordinates and the following 
*	    formula.  distance = sqrt(((x1-x2)^2) + ((y1-y2)^2))
*/

#include <stdio.h>
#include <math.h>

/* function prototypes */
double distance(double x1, double x2, double y1, double y2);

int main()
{
	double length;
	double x1, x2, y1, y2;
	
	printf("\n\n");
	printf("You will now enter the cartesian coordinates for two points");
	printf("\nPlease enter values for x1 and y1 > ");
	scanf("%lf%lf", &x1, &y1);
	printf("Please enter values for x2 and y2 > "); 
	scanf("%lf%lf", &x2, &y2);

	length = distance(x1, x2, y1, y2);
	printf("The distance between these points is: %.1f\n\n", length);

	return(0);
}
/* Function to calculate the distance between points based on their
*  cartesian coordinates */

double
distance(double x1, double x2, double y1, double y2)
{
	double distance;
	
	distance = sqrt((pow(x1-x2, 2)) + (pow(y1-y2, 2)));
	return(distance);
}
