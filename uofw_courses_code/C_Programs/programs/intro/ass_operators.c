/* ass_operators */
/* Author: Malachi Griffith 
*  Date: Sept. 17 2002
*  Purpose: Illustrate use of assignment operators
*/

/* The following operators are defined as follows:
*	+=  means "is increased by"
*	-=  means "is decreased by"
*	(asterix)=  means "is multiplied by"
*	(slash)=  means "is divided by"
*/

/* example: A program to calculate the balance of a bank account */
/* Note that assignment operators have considerably lower precedence than
*  the arithmetic operators. */

#include <stdio.h>

main()
{
	double balance;
	double deposits;
	double withdrawals;

	printf("Please type in the opening balance:\n");
	scanf("%lf", &balance);
	
	printf("Please type in the deposit amount:\n");
	scanf("%lf", &deposits);

	balance += deposits;

	printf("Please type in the withdrawal amount:\n");
	scanf("%lf", &withdrawals);

	balance -= withdrawals;

	printf("The final balance in the account is $%.2f\n", balance);
}
