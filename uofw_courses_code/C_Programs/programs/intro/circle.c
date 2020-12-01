/* 
* Purpose: Calculate the area and circumferance of a circle
* Author: Malachi Griffith
* Date: Sept. 13 2002
*/

#include <math.h>
#include <stdio.h>
#define pi 3.14159

int
main(void)
{
	double radius,
	       area,
               circumf;

	int    num_circ;

	char   circ_name;


	printf("Enter a single letter name for the circle>");
	scanf("%c", &circ_name);
	printf("Enter the radius of this circle>");
	scanf("%lf", &radius);

	area = pi * (radius * radius);
	circumf = 2 * pi * radius;

	printf("Area of circle %c with radius %f is %f. \n",
		circ_name,radius, area);
	printf("Circumferance of circle %c with radius %f is %f. \n",
		circ_name, radius, circumf);

return(0);
}
