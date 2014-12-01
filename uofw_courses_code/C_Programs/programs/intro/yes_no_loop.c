/*yes_no_loop.c*/
/*
* Author: Malachi Griffith
* Date: Sept. 15 2002
* Purpose: Calculates the dollar value for various combinations of 
*          loose change entered by the user.
*/

#include <stdio.h>

main()

{
	/* Program Inputs */
	char first, middle, last;  /* Customer's initials */
	int quarters;  /* Number of quarters */
	int dimes;  /* Number of dimes */
	int nickels;  /* Number of nickels */
	int pennies;  /* Number of pennies */

	/* Program Outputs */
	int dollars;  /* Equivalent number of dollars */
	int change;  /* Change left over */

	/* Other Program Variables */
	int total_cents;  /* Total number of cents */
	char repeat = 1; /* Variable to allow repeat of program */
	char answer;
	/* Collect the input data from the user */

while (repeat == 1)
	{
	
	printf("\nPlease enter your first, middle, and last initials >");
	scanf(" %c %c %c", &first, &middle, &last);

	printf("How many quarters do you have? >");
	scanf(" %d", &quarters);
	
	printf("How many dimes do you have? >");
	scanf(" %d", &dimes);

	printf("How many nickels do you have? >");
	scanf(" %d", &nickels);
	
	printf("How many pennies do you have ? >");
	scanf(" %d", &pennies);             

	/* Calculate the total cents for the user */

	total_cents = (quarters * 25) + (dimes * 10) + 
                      (nickels * 5) + (pennies);

	/* Calculate the value of this in dollars and cents */

	dollars = total_cents / 100;
	change = total_cents % 100;

	/* Display the result */
	
	printf
	("\n\nThe value of your coin collection is %d dollars and %d cents\n",
	dollars, change);

	printf("Would you like to calculate another collection? (y or n) >");
	scanf(" %c", &answer);
	if (answer == 'y') repeat = 1;
	else repeat = 0;
	}
}

