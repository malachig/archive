/* quadrant.c */

/* Author: Malachi Griffith
*  Date: Oct. 25 2002
*  Purpose: Based on a set of cartesian coordinates this program tells
*  you which quadrant they lie in;
*		QII  |  QI
*		----------
*	  	QIII |  QIV
*/

#include <stdio.h>

main()
{
	double x, y;

	printf("\nEnter the x and y components seperated by a space > ");
	scanf("%lf%lf", &x, &y);

	if (x == 0)
		printf("\nThe point is on the y axis\n\n");
	else if (y == 0)
		printf("\nThe point is on the x axis\n\n");

	else if (x > 0 && y > 0)
		printf("\nThe point is in quadrant I\n\n");
	else if (x > 0 && y < 0)
		printf("\nThe point is in quadrant IV\n\n");
	else if (x < 0 && y > 0)
		printf("\nThe point is in quadrant II\n\n");
	else if (x < 0 && y < 0)
		printf("\nThe point is in quadrant III\n\n");
}
