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
        printf("debugging  %f  \n", loan_amount);
	fscanf(ifptr, "%lf", &annual_income);
        printf("debugging  %f  \n", annual_income);

	{
	number++;
	if ((loan_amount < 5000.0) && (annual_income >= 30000.0))
		printf("\n\nLoan is approved!\n");

	fscanf(ifptr, "%f", &loan_amount);
	fscanf(ifptr, "%f", &annual_income);
	}
	fclose(ifptr);
	return(0);
}



