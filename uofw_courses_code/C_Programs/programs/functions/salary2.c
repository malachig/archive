/* Author: Malachi Griffith
*  Date: Oct. 8 2002
*  Purpose: Computes the users gross salary based on wage and hours.
*/

#include <stdio.h>

int
main(void)
{
	double 	wage,
		hours,
		overtime,
		bonus,
		salary;

	printf("\nPlease enter your hourly wage > ");
	scanf("%lf", &wage);
	printf("Please enter the hours worked > ");
	scanf("%lf", &hours);
	printf("Please enter the number of overtime hours > ");
	scanf("%lf", &overtime);

	salary = wage * hours;
	bonus = overtime * (wage * 1.5);
	salary += bonus;	

	printf("The salary resulting is: %.2f\n\n", salary);

	return(0);
}



