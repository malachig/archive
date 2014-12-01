/* dayofyear.c */

/* Author: Malachi Griffith
*  Date: Oct. 25 2002 
*  Purpose: Takes month, day and year and calculates the day of the year 
*  (1 to 366) 
*/

#include <stdio.h>

/* Function prototype */
int leap(int year_value);

main()
{
	int day,
	    month,
	    year;

	int days_upto_month[13] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
				   304, 334, 365};
	int day_of_year;
	int extra_day;

	printf("\nEnter the month, day and year, seperated by spaces > ");
	scanf("%d%d%d", &month, &day, &year);

	day_of_year = (days_upto_month[month - 1]) + (day);	

	extra_day = leap(year);	
	
	if (month > 2)
		day_of_year += extra_day;

	printf("\nThat day of the year is %d\n\n", day_of_year);
}

int
leap(int year_value)
{
	int remainder4;
	int remainder100;
	int remainder400;
	int extra;

	remainder4 = year_value % 4;
	remainder100 = year_value % 100;
	remainder400 = year_value % 400;	

	if (remainder4 == 0 && remainder100 != 0)
		extra = 1;
	else if (remainder4 == 0 && remainder100 == 0 && remainder400 == 0)
		extra = 1;
	else
		extra = 0;

	return (extra);
}
		
