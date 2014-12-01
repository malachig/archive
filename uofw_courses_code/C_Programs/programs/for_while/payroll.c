/* payroll.c */

/* Author: Malachi Griffith
*  Date:  Oct. 5 2002
*  Purpose: Compute the payroll of a company 
*/

#include <stdio.h>

int
main()
{
	double total_pay;		/* company payroll */
	int count_emp;			/* current employee */
	int number_emp;			/* number of employees */
	double hours;			/* hours worked */
	double rate;			/* hourly rate */
	double pay;			/* pay for this period */

	/* Get the number of employees. */
	printf("Enter number of employess. > ");
	scanf("%d", &number_emp);

	/* Compute each employee's pay and add it to the payroll. */
	total_pay = 0.0;
	count_emp = 0;
	
	while (count_emp < number_emp)
	{
		printf("Hours > ");
		scanf("%lf", &hours);
		printf("Rate > $");
		scanf("%lf", &rate);
		pay = hours * rate;
		printf("Pay is %6.2f\n\n", pay);
		total_pay = total_pay + pay;
		count_emp++;
	}
	
	printf("All employees processed\n");
	printf("Total payroll is $%8.2f\n", total_pay);

	return(0);
}



