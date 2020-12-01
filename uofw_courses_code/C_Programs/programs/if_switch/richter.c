/* richter.c */

/* Author: Malachi Griffith
*  Date: Oct. 25 2002
*  Purpose: Determines damage depending on strength of earthquake.
*/

#include <stdio.h>

main()
{
	double richter;
	
	printf("\nPlease enter the value on the richter scale > ");
	scanf("%lf", &richter);

	if (richter < 5.0)
		printf("\nLittle or no damage\n\n");
	else if(richter < 5.5)
		printf("\nSome damage\n\n");
	else if (richter < 6.5)
		printf("\nSerious damage, walls may crack or fall\n\n");
	else if (richter < 7.5)
		printf("\nDisaster: houses and buildings may collapse\n\n");
	else if (richter >= 7.5)
		printf("\nCatastrophe: pray for lives\n\n");
}

