/* airplane.c */

/* Author: Malachi Griffith
*  Date: Oct 2 2002
*  Purpose: To assign a classification to planes based on their speed and size.
*/

#include <stdio.h>

main()
{
	int speed;
	int length;

	printf("\n\nPlease enter the speed in Km/Hr as an integer > ");
	scanf(" %d", &speed);
	printf("Please enter the length of the plane in meters as an integer> ");
	scanf(" %d", &length);

	if (speed > 1100)
		if (length > 52)
			printf("\nThe plane is civilian in nature.\n\n");
		else 
			printf("\nThe plane is military in nature.\n\n");
	else
		printf("\nThe aircraft type is unknown.\n\n");
	return(0);
}



