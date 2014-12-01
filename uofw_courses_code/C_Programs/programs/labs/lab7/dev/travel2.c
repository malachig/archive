/* Author: Malachi Griffith
*  Date:  Oct. 16 2002
*  Purpose: Use an array to store ticket prices for a flight to one
*	    of 6 different cities.  Get the city from the use rand the 
*	    number of tickets they want and calculate the cost.  Display
*	    an error message for invalid city number.
*/

#include <stdio.h>
#define MAX_CITIES 6

/* Function prototypes */
void total_tickets(int num_tickets, int city, int *tickets);

main()
{
	double cost[] = {56.79, 104, 93.4, 14.6, 84.8, 74.9};
	int tickets[MAX_CITIES] = {0};

	int city;
	int num_tickets = 0;
	double cost_for_city = 0;
	double total_cost = 0;

	FILE *ifptr;

	ifptr = fopen ("travel.dat", "r");

	fscanf(ifptr, "%d", &city);
	fscanf(ifptr, "%d", &num_tickets);

	while(!feof(ifptr))
	{
	total_tickets(num_tickets, city, &tickets[city]);
	
	fscanf(ifptr, "%d", &city);
	fscanf(ifptr, "%d", &num_tickets);
	}			

	for (city = 1; city <= 6; city++)
	{
		cost_for_city = tickets[city - 1] * cost[city - 1];
		total_cost += cost_for_city;	
		printf("\nThe total tickets for city %d is %d",
			 city, tickets[city - 1]);
		printf("\nThe total cost for these tickets is $%.2f", 
			cost_for_city);	
	}

		printf("\nThe total cost of your purchase is $%.2f\n\n",
	         	 total_cost);

	fclose(ifptr);
	return(0);
}

/* Function total_tickets */
void 
total_tickets(int num_tickets, int city, int tickets[])
{
	tickets[city - 1] += num_tickets;
}





