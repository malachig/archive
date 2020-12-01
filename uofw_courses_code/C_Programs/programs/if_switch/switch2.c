/* switch2.c */

/* Author: Malachi Griffith
*  Date: Oct. 2 2002 
*  Purpose: Use the switch statement to assign a brightness value
*	    to lightbulbs of varying powers.
*/

#include <stdio.h>

main()
{
	int lumens = 0;
	int watts;

	printf("\n\nPlease enter the Watts of your bulb (15,25,40,60,75,100) > ");
	scanf(" %d", &watts);

	switch (watts)
	{
	case 15:
		lumens = 125;
		break;
	
	case 25:
		lumens = 215;
		break;
	
	case 40:
		lumens = 500;
		break;
	
	case 60:
		lumens = 880;
		break;
	
	case 75:
		lumens = 1000;
		break;
	
	case 100:
		lumens = 1675;
		break;

	default:
		printf("\nNot a valid bulb type!\n\n");
	}
	if (lumens != 0)
	printf("\nThe lumens associated with a bulb of %d watts is: %d\n\n", watts, lumens);





}



