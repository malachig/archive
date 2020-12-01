/* Author: Malachi Griffith
*  Date:  Oct. 16 2002
*  Purpose: Use an array to store ticket prices for a flight to one
*	    of 6 different cities.  Get the city from the use rand the 
*	    number of tickets they want and calculate the cost.  Display
	    an error message for invalid city number.
*/

#include <stdio.h>

main()
{
	double city[] = {56.79, 104, 93.4, 14.6, 84.8, 74.9};

	int i;

	int num_tickets;
	
	double total_cost;

	printf("\nPlease enter the number of tickets you want > ");
	scanf("%d", &num_tickets);

	printf("Which city would you like to fly to (1 to 6)? > ");
	scanf("%d", &i);
	
	if (i < 1 || i > 6)
	{
		printf("\nInvalid city entry!");
		exit();
	}
	else
	{
		total_cost = (double) (num_tickets * city[i - 1]);

		printf("\nThe total cost of your purchase is $%.2f\n\n",
	         	 total_cost);
	}
}
