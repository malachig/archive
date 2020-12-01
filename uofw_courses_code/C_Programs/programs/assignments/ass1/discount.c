/* Assignment 1: Question 3 Part A. (7 Marks) discount.c */

/* Author: Malachi Griffith
*  Date: Sept. 28 2002
*  Purpose: To ask the user for the price of an order, determine the 
*	    appropriate discount using if/else statements, and calculate 
*	    the net cost of the order.
*/

#include <stdio.h>

main()
{
	/* define the variables */	

	double price,  		/* Price of order */
	       discount_rate,	/* Rate of discount, depending on price */
	       discount,	/* Price multiplied by discount */
	       net_price;	/* Price when discount is taken into account */

	/* Retrieve price from the user */	

	printf("\n\nPlease enter the total price of your order > ");
	scanf(" %lf", &price);


	/* Determine which discount rate applies based on price.
	*  When the price is less than $200.00, the discount is only,
	*  10%.  If the price is greater than $200.00, a discount of 
	*  15% applies. 
	*/
	
	if (price <= 200.0)
		discount_rate = 0.10;
	else
		discount_rate = 0.15;


	/* Calculate the net price and display results */

	discount = price * discount_rate;
	net_price = price - discount;

	printf("\nThe price of the order was: $%.2f", price);
	printf("\nThe discount amount was: $%.2f (at %.0f percent)",
	       discount, discount_rate * 100);
	printf("\nThe net price after the discount is: $%.2f\n\n",
	        net_price);
	
	return(0); 
}



