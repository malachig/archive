/* payroll5.c */

/* Author: Malachi Griffith
*  Date:  Oct. 5 2002
*  Purpose: Compute the payroll of a company 
*/

#include <stdio.h>

int
main()
{
	double total_pay;		/* company payroll */
	double hours;			/* hours worked */
	double rate;			/* hourly rate */
	double pay;			/* pay for this period */

	FILE * ifptr;
	ifptr = fopen("payroll.dat", "r");


	/* Compute each employee's pay and add it to the payroll. */
	total_pay = 0.0;
		fscanf(ifptr, "%lf", &hours);
		fscanf(ifptr, "%lf", &rate);
	
	while (!feof(ifptr))
	{
		printf("Hours are: %.0f\t", hours);
		printf("Rate is: $%.2f\t", rate);
		pay = hours * rate;
		printf("Pay is %6.2f\n", pay);
		total_pay = total_pay + pay;
	
		fscanf(ifptr, "%lf", &hours);
		fscanf(ifptr, "%lf", &rate);
	}
	
	printf("All employees processed\n");
	printf("Total payroll is $%8.2f\n", total_pay);

	fclose(ifptr);
	return(0);
}



