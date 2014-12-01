/* payroll4.c */

/* Author: Malachi Griffith
*  Date:  Oct. 5 2002
*  Purpose: Compute the payroll of a company 
*/

#include <stdio.h>
#define SENTINEL -99

int
main()
{
	double total_pay;		/* company payroll */
	int count_emp;			/* current employee */
	double pay;			/* pay for this period */

	/* Compute each employee's pay and add it to the payroll. */
	total_pay = 0.0;
	printf("\n\nEnter -99 to quit entry");
	printf("\nEnter the first pay value > $");
	
	for (scanf("%lf", &pay);	
	     pay != -99;
	     scanf("%lf", &pay))	
	{
		printf("Pay is %6.2f\n", pay);
		printf("Pay > $");
		total_pay = total_pay + pay;
	}
	
	printf("All employees processed\n");
	printf("Total payroll is $%8.2f\n", total_pay);

	return(0);
}



