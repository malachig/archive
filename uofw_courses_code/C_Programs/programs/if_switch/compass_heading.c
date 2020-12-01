/* compass_heading.c */

/* Author: Malachi Griffith
*  Date: Oct. 24 2002
*  Purpose: Convert a degree value between 0 and 360 to a compass heading.
*/

#include <stdio.h>

main()
{
	double degrees;
	double degrees_to_turn;	
	char face; /* N for north, S for South */
	char turn; /* E for East, W for West */


	printf("Enter the degree value > ");
	scanf("%lf", &degrees);
	while (degrees > 360 || degrees < 0)
	{
		printf("Not a valid entry, try again!");	
		printf("\nEnter the degree value > ");
		scanf("%lf", &degrees);
	}

	if (degrees >= 270 && degrees <= 360)
	{
		face = 'N';
		turn = 'W';
		degrees_to_turn = 360 - degrees;
	}
	else if (degrees >= 0 && degrees <= 90)
	{
		face = 'N';
		turn = 'E';
		degrees_to_turn = degrees;
	}

	else if (degrees > 90 && degrees <= 180)
	{
		face = 'S';
		turn = 'E';
		degrees_to_turn = 180 - degrees;
	}

	else if (degrees > 180 && degrees < 270)
	{
		face = 'S';
		turn = 'W';
		degrees_to_turn = degrees - 180;
	}

	if (face == 'S')
		printf("\nFace South");
	else 
		printf("\nFace North");
	
	printf(", then turn %.0f degrees", degrees_to_turn);

	if (turn == 'W')
		printf(", to the West\n\n");
	else
		printf(", to the East\n\n");
}
