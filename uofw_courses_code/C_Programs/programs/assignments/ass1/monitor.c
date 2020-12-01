/* Assignment 1: Question 1 (10 Marks)  monitor.c */

/* Author: Malachi Griffith
*  Date: Sept. 23 2002 
*  Purpose: Receives orders for monitors of three different types from 
*           the user and calculates the cost of the monitors by type,  
*	    as well as finding the total cost of the order.
*/

#include <stdio.h>
main()
{
	/* Define the Variables */
	double 	cost_a = 500.0,		/* Cost of Monitor Type A */
		cost_b = 695.5,		/* Cost of Monitor Type B */
		cost_c = 799.99;	/* Cost of Monitor Type C */

	int	monitors_a,		/* Number of A Monitors ordered */
		monitors_b,		/* Number of B Monitors ordered */
		monitors_c;  		/* Number of C Monitors ordered */

	double 	total_a,		/* Total Cost of A Type Monitors */
		total_b,		/* Total Cost of B Type Monitors */
		total_c,		/* Total Cost of C Type Monitors */
		grand_total;		/* Total cost of the entire order */

	/* Get Input from User - Number of each monitor type wanted */
	printf("\n\n");
	printf("How many Type A Monitors would you like to order? > ");
	scanf("%d", &monitors_a);

	printf("How many Type B Monitors would you like to order? > ");
	scanf("%d", &monitors_b);

	printf("How many Type C Monitors would you like to order? > ");	
	scanf("%d", &monitors_c);

	/* Echo the order back to the user for confirmation */	
	printf
	("\nOrder is: %d Type A, %d Type B, and %d Type C Monitors.\n",
        monitors_a, monitors_b, monitors_c);

	/* Calculate the cost for each type of montitors and diplay it */
	total_a = monitors_a * cost_a;
	printf("\nCost for the Type A Monitors is $%.2f", total_a);

	total_b = monitors_b * cost_b;
	printf("\nCost for the Type B Monitors is $%.2f", total_b);

	total_c = monitors_c * cost_c;
	printf("\nCost for the Type C Monitors is $%.2f\n", total_c);

	grand_total = total_a + total_b + total_c;
	printf("\nThe grand total of your order is: $%.2f\n\n",
	       grand_total); 	
	
	return(0);
}



