/* Author: Malachi Griffith
*  Date: Oct. 12 2002
*  Purpose: Reads account info from a data file called bank.dat and
*	    unpdates a bak account.  A value can either be W for withdrawal
*	    D for deposit, or Q for quit (sentinel value).  The first number
*	    in the datafile is the initial balance.
*/

#include <stdio.h>

main()
{
	char action;

	double value;
	double balance;

	FILE * input_file;

	input_file = fopen("bank.dat", "r");

	/*Read in the initial balance */
	fscanf(input_file, "%lf", &balance);
	printf("\nThe initial account balance is: $%.2f", balance);

	fscanf(input_file, " %c", &action);
	fscanf(input_file, " %lf", &value);

	while (action != 'Q')
	{
	if (action == 'W')
		balance -= value;
	else if (action == 'D')
		balance += value;
	else
		printf("\nInvalid entry in datafile.\n");
	fscanf(input_file, " %c", &action);
	fscanf(input_file, " %lf", &value);
	}
	printf("\nThe final account balance is $%.2f\n", balance);	
	
	fclose(input_file);
	return(0);
}



