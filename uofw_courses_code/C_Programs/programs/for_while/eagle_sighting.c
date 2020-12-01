/* eagle_sighting.c */

/* Author: Malachi Griffith
*  Date: Oct. 5 2002
*  Purpose: Tally by month the bald eagle sightings for the year.
*	    Each month's sightings are terminated by the sentinel
*	    zero.
*	NOTE:  To enter data into this program enter the sightings 
*	of each member at the prompt, seperated by spaces.  When all the 
* 	sightings for a particular month are entered, enter the sentinel
* 	value of 0 to continue to the next month.
*/

#include <stdio.h>

#define SENTINEL 0
#define NUM_MONTHS 12

int
main(void)
{
	int month,	/* number of month being processed */
	mem_sight,	/* one member's sightings for this month */
	sightings;	/* total sightings so far for this month */

	printf("BALD EAGLE SIGHTINGS\n");
	for (month = 1;
             month <= NUM_MONTHS;
	     ++month)
	{
		sightings = 0;
		scanf("%d", &mem_sight);
	
		while (mem_sight != SENTINEL)
		{
			if (mem_sight >= 0)
				sightings += mem_sight;
			else 
				printf("Warning, negative count %d ignored\n",
					mem_sight);
			scanf("%d", &mem_sight);
		} /* inner while */
	
	printf("  month %2d: %2d\n", month, sightings);
	} /* outer for */

	return(0);
}



