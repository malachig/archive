/* Author: Malachi Griffith
*  Date: Oct 14 2002
*  Purpose: A program to calculate shipping charges by region. 
*/

#include <stdio.h>
#define MAX_RATE 5

main()
{
	int region;
	
	double rate[MAX_RATE] = {0.075, 0.080, 0.082, 0.085, 0.088};

	double price;
	double charge;
	double total_price;

	printf("Enter the price of the item: $");
	scanf("%lf", &price);
	printf("\nEnter shipping region (1-5): ");
	scanf("%d", &region);
	printf("\n");

	if ((region < 1) || (region > 5))
		printf("Invalid region number entered.  Try again.\n");
	else
	{
		charge = price * rate[region - 1];
		total_price = price + charge;
		printf("\n\nItem Price:   %9.2f", price);
		printf("\nShipping Charge: %9.2f", charge);
		printf("\nTotal Price:   %9.2f", total_price);
		printf("\n");
	}
}



