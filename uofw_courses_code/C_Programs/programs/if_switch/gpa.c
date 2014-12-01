/* gpa.c */

/* Author: Malachi Griffith
*  Date:  Oct. 2 2002
*  Purpose: Use of multiple alternative decisions with if-else
*/

#include <stdio.h>

main()
{
	double gpa;

	printf("\n\nPlease enter your GPA for the semester > ");	
	scanf(" %lf", &gpa);
	
	if (gpa < 0.0 || gpa > 4.00)
		printf("\n\nNot a valid GPA!\n");
	else if(gpa <= 0.99)
		printf("\nFailed semester - registration suspended\n");
	else if (gpa <= 1.99)
		printf("\nOn probation for next semester\n");
	else if (gpa <= 2.99)
		printf("\nDoin fine baby ...\n");	
	else if (gpa <= 3.49)
		printf("\nDean's list for that semester\n");
	else 
		printf("\nHighest honours for that semester\n");			
}
