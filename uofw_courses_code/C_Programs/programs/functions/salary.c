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
		salary;

	printf("\nPlease enter your hourly wage > ");
	scanf("%lf", &wage);
	printf("\nPlease enter the hours worked > ");
	scanf("%lf", &hours);

	salary = wage * hours;
	
	printf("The salary resulting is: %.2f\n\n", salary);

	return(0);
}



