/* Author: Malachi Griffith
*  Date: Nov 6 2002 
*  Purpose: A program to demonstrate gets and string output.  Based on 
*  Molluzo, p 308.
*  The gets() function retrives a string from the standard input and places
*  into an array already defined in the program.  Make sure this array is
*  long enough to fit the entire string.
*/

#include <stdio.h>

#define TAX_RATE 0.07
#define ARRAY_SIZE 51

main()
{
	char item_name[ARRAY_SIZE];

	double price;
	double tax;
	double total;

	printf("\nEnter the name of the item: ");
	gets(item_name);

	printf("\nEnter item price: ");
	scanf("%lf", &price);

	tax = price * TAX_RATE;
	total = price + tax;

	printf("\nItem: %s\n", item_name);
	printf("\nPrice: %9.2f", price);
	printf("\nTax:  %9.2f", tax);
	printf("\n%s", "----------------");
	printf("\n%-7s%9.2f\n", "Total:", total);
	printf("\n\n");
}
