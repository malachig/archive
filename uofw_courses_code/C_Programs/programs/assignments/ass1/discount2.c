/* Assignment 1: Question 3 Part B. (8 Marks) discount2.c */

/* Author: Malachi Griffith
*  Date: Sept. 28 2002
*  Purpose: To ask the user for the price of an order, determine the 
*	    appropriate discount using if/else statements, and calculate
*	    the net cost of the order.
*/

#include <stdio.h>

main()
{
	/*  Define the variables */

	double price,  		/* Price of order */
	       discount_rate,	/* Rate of discount, depending on price */
	       discount,	/* Price multiplied by discount */
	       net_price;	/* Price when discount is taken into account */
	int age;

	/* Retrieve price from the user */	

	printf("\n\nPlease enter the total price of your order > ");
	scanf(" %lf", &price);
	printf("Please enter your age > ");
	scanf(" %d", &age);


	/* Determine which discount rate applies based on price.
	*  When the price is less than or equal to $200, seniors 
	*  (65 or older) get a discount of 15% and others get %10.
	*  When the price is greater than $200, seniors get a discount
	*  of 20% and everyone else gets 15%.
	*/
	
	if (price <= 200.0)
		if (age >= 65)	
			discount_rate = 0.15;
		else
			discount_rate = 0.10;
	
	else
		if (age >= 65)
			discount_rate = 0.20;
		else
			discount_rate = 0.15;

	/* Calculate the net price and display results */
	
	discount = price * discount_rate;
	net_price = price - discount;

	printf("\nThe price of the order was: $%.2f", price);
	printf("\nYour discount rate is: %3.0f percent", discount_rate * 100);
	printf("\nThe discount amount was: $%.2f", discount);
	printf("\nThe net price after the discount is: $%.2f\n\n", net_price);
	
	return(0); 
}



