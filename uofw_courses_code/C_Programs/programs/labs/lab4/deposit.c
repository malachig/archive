/*  Author: Malachi Griffith
*   Date: Sept. 26 2002
*   Purpose: Use a loop to calculate how many months it will be before 
*	     an account exceeds $1000 when Each month the deposits follow
*	     the pattern: $1, $2, $4, etc.
*/

#include <stdio.h>
#include <math.h>
main()

{
	double
	       increase,
 	       balance,
	       month;
	month = 1.0;
	balance = 1.0;
	printf("\n\nBalance after one month is: $%.2f\n", balance);
	while (balance < 1000.0)
	{	
	increase = pow(2,month);
	balance += increase;  	
	month++;
	printf("Balance after month %.0f is $%.2f\n", month, balance);	
	}

	printf("Number of months it took for balance to exceed $1000, is %.0f
	\n\n",month);
}








