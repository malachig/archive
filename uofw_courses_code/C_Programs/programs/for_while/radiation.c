/* radiation.c */

/* Author: Malachi Griffith
*  Date: Oct. 6 2002
*  Purpose: Calculates and displays a table showing the safety level 
*	    of a coffee room.
*/

#include <stdio.h>

#define SAFE_RAD 0.466	/* Safe level of radiation */
#define SAFETY_FACT 10.0	/* Safety factor */

int rad_table(double init_radiation, double min_radiation);

int
main(void)
{
	int day;	/* day user can enter room */
	double init_radiation; /* radiation level right after leak */
	double min_radiation;	/*safe level divided by safety factor */

	/* Compute stopping level of radiation. */
	min_radiation = SAFE_RAD / SAFETY_FACT;

	/* Prompts user to enter initial radiation level */
	printf("Enter the radiation level (in millirems) > ");
	scanf("%lf", &init_radiation);

	/* Displays table */
	day = rad_table(init_radiation, min_radiation);

	/* Displays day the user can enter the room */
	printf("\nYou can enter the room on day %d. \n", day);

	return(0);
}

/* Diplays a table showing the radiation level and safety status every 3
* days until the room is deemed safe to enter. Reurns the day number for 
* the fist safe day. 
* Pre: min_radiation and init_radiation are defined.
* Post: radiation_lev <= min_radiation
*/

int
rad_table(double init_radiation, double min_radiation)
{
	int day;  /* days elapsed since substance leak */
	double radiation_lev;	/* current radiation level */

	day = 0;

	printf("\n   Day   Radiation   Status\n      (millirems)\n");
	for (radiation_lev = init_radiation;
	     radiation_lev > min_radiation;
	     radiation_lev /= 2.0)
	{
	if (radiation_lev > SAFE_RAD)
		printf("  %3d%5c%9.4f   Unsafe\n", day, ' ', radiation_lev);
	else
		printf("  %3d%5c%9.4f   Safe\n", day, ' ', radiation_lev);
	day += 3;
	}

	return(day);
}



