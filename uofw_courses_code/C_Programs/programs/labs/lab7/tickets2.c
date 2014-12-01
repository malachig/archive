/* Author: Malachi Griffith
*  Date:  Oct. 16 2002
*  Purpose: Use an array to store ticket prices for a flight to one
*	    of 6 different cities.  Get the city from the use rand the 
*	    number of tickets they want and calculate the cost. 
*/

#include <stdio.h>
#define MAX_CITIES 6

/* Function prototypes */
void total_tickets(int num_tickets, int city, int *tickets);
void ticket_cost(int city, int tickets[], double cost[],
	         double cost_per_city[], double *total_cost);
void print_out(int tickets[], double cost_per_city[], double total_cost); 

main()
{
	/* Array that defines the cost of a flight to each city */
	double cost[] = {56.79, 104, 93.4, 14.6, 84.8, 74.9};

	/* Array to tally tickets purchased for each city */	
	int tickets[MAX_CITIES] = {0};

	/* Array to tally the total cost of tickets to each city */
	double cost_per_city[MAX_CITIES] = {0};

	int city;
	int num_tickets = 0;
	double total_cost;

	total_tickets(num_tickets, city, tickets);
	ticket_cost(city, tickets, cost, cost_per_city,
		    &total_cost);
	print_out(tickets, cost_per_city, total_cost);

	return(0);
}

/*
*  Function:  total_tickets 
*  Adds up the tickets purchased for each city.
*/
void 
total_tickets(int num_tickets, int city, int tickets[])
{
	FILE *ifptr;
	ifptr = fopen ("travel.dat", "r");
	
	fscanf(ifptr, "%d", &city);
	fscanf(ifptr, "%d", &num_tickets);
	
	while(!feof(ifptr))
	{
	tickets[city - 1] += num_tickets;
	
	fscanf(ifptr, "%d", &city);
	fscanf(ifptr, "%d", &num_tickets);
	}
	
	fclose(ifptr);
}

/*
*  Function: ticket_cost
*  Calculates the cost of the tickets for each city.  
*  Pre: The cost and tickets arrays are already defined.
*  Post: It updates the values of the cost_per_city array.
*/
void 
ticket_cost(int city, int tickets[], double cost[], double cost_per_city[],
	    double *total_cost)
{
	double total;
	
	for (city = 1; city <= 6; city++)
		{
		cost_per_city[city-1] = tickets[city - 1] * cost[city - 1];
		total += cost_per_city[city -1];	
		}
	*total_cost = total;
}

/*
*  Function: print_out
*  Prints out the resulting cost per city, the total tickets to each city 
*  and the the grand total cost of the order.
*  Pre:  The array tickets and cost_per_city are defined.
*/ 

void 
print_out(int tickets[], double cost_per_city[], double total_cost) 
{
	int city;	
	for (city = 1; city <= 6; city++)
	{
		printf("\nThe total tickets for city %d is %d",
			 city, tickets[city - 1]);
		printf("\nThe total cost for these tickets is $%.2f", 
			cost_per_city[city -1]);	
	}

		printf("\nThe total cost of your purchase is $%.2f\n\n",
	         	 total_cost);
}
