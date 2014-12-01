/* Author: Malachi Griffith
*  Date: Oct. 3 2002
*  Purpose: To screen loan applications from a datafile
*/

#include <stdio.h>

main()
{

	double loan_amount;
	double annual_income;
	int number = 0;

	FILE * ifptr;
	
	ifptr = fopen("application.dat", "r");

	fscanf(ifptr, "%lf", &loan_amount);
	fscanf(ifptr, "%lf", &annual_income);

	while(!feof(ifptr))
	{
	number++;
	if (loan_amount < 5000.0 && annual_income >= 30000.0)
		printf("\n\nLoan is approved!\n");
	else if (loan_amount >= 5000.0 && loan_amount < 20000.0 
		 && annual_income >= 75000.0)
		printf("\n\nLoan is approved!\n");
	else
		printf("\n\nLoan is rejected!\n");

	fscanf(ifptr, "%lf", &loan_amount);
	fscanf(ifptr, "%lf", &annual_income);
	}
	fclose(ifptr);
	return(0);
}



