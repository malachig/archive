/* payroll3.c */

/* Author: Malachi Griffith
*  Date:  Oct. 5 2002
*  Purpose: Compute the payroll of a company in a subfunction and
*	    return the total payroll value to the main function.
*	    The main function gets the initial input from the user
*	    and passes on the number of employees to the function,
*	    'payroll'.
*/

#include <stdio.h>

/* Function prototype */
double payroll(int number_emp);


int
main(void)
{
	double total;
	int number_emp;			/* number of employees */

	/* Get the number of employees. */
	printf("Enter number of employess. > ");
	scanf("%d", &number_emp);

	total = payroll(number_emp);
	
	printf("All employees processed\n");

	printf("Total payroll is $%8.2f\n", total);
	return(0);
}


/*
* FUNCTION Compute each employee's pay and add it to the payroll.
*/
double
payroll(number_emp)
{
	double total_pay;		/* company payroll */
	int count_emp;			/* current employee */
	double hours;			/* hours worked */
	double rate;			/* hourly rate */
	double pay;			/* pay for this period */

	total_pay = 0.0;
	count_emp = 0;
	
	for (count_emp = 0;		/* initialisation	*/
	     count_emp < number_emp;	/* repetition condition	*/
	     count_emp += 1)		/* update		*/
	{
		printf("Hours > ");
		scanf("%lf", &hours);
		printf("Rate > $");
		scanf("%lf", &rate);
		pay = hours * rate;
		printf("Pay is %6.2f\n\n", pay);
		total_pay = total_pay + pay;
	}
	return (total_pay);
}	




