/*yard_area.c*/
/* Author: Malachi Griffith
*  Date: Sept. 29 2002 
*  Purpose: Calculate the time required to mow a yard based on the size of the yard and 
*	    a house situated in the yard.
*/

#include <stdio.h>

main()
{
	int mow_rate = 2,  /* Grass can be moved at 2 sq. ft. per second */
	    yard_width,
	    yard_length, 
	    yard_area,
	    house_width,
	    house_length, 
	    house_area,
	    mow_area;
	double mow_time;

	printf("\n\nEnter the width of the yard in feet > ");
	scanf(" %d", &yard_width);
	printf("Enter the length of the yard in feet > ");
	scanf(" %d", &yard_length);
	printf("Enter the width of the house in feet > ");
	scanf(" %d", &house_width);
	printf("Enter the length of the house in feet > ");
	scanf(" %d", &house_length);

	yard_area = yard_width * yard_length;
	house_area = house_width * house_length;	
	mow_area = yard_area - house_area;
	mow_time = (double) (mow_area / mow_rate);

	printf("The time required to mow this are is: %7.2lf seconds\n\n", mow_time);

}
