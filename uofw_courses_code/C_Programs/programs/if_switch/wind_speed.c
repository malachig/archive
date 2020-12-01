/* wind_speed.c */

/* Author: Malachi Griffith
*  Date: Sept. 2 2002
*  Purpose: Use multiple-alternative if statement to code the following
*/

#include <stdio.h>

main()
{
	int wind;

	printf("\n\nPlease enter the wind speed as an integer > ");
	scanf(" %d", &wind);

	if (wind < 0)
		printf("\nNot a valid wind speed\n\n");
	else if (wind < 25)
		printf("\nNot a strong wind\n\n");
	else if (wind <= 38)
		printf("\nStrong wind\n\n");
	else if (wind <=54)
		printf("\nGale\n\n");
	else if (wind <= 72)
		printf("\nWhole gale\n\n");
	else 
		printf("\nHurricane\n\n");
return(0);
}



