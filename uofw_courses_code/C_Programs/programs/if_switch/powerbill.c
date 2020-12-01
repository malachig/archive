/* powerbill.c */

/* Author: Malachi Griffith
*  Date: Oct. 25 2002
*  Purpose: Calculates power bill depending on the type of user and 
*  their usage.
*  R = residential user, C = commercial user, and I = industrial user.
*/

#include <stdio.h>

main()
{
	char user_type;
	int account_number;
	int power_used;
	double total_bill;
	double extra_charge;

	printf("\nEnter your account number > ");
	scanf("%d", &account_number);
	printf("Enter your code (R, C, or I) > ");	
	scanf(" %c", &user_type);
	printf("Enter the kilowatt hours of power used as whole number > ");
	scanf("%d", &power_used);

	if (user_type == 'R' || user_type == 'r')
		total_bill = (double) (6.00 + (0.052 * power_used));
	
	else if (user_type == 'C' || user_type == 'c')
		{
			total_bill = 60.00;
			if (power_used > 1000)
			{
				extra_charge = (power_used - 1000) * 0.045;
				total_bill += extra_charge;
			}
		}

	else if (user_type == 'I' || user_type == 'i')
		{
			total_bill = 76.00;
			if (power_used > 1000)
			{
				extra_charge = (power_used - 1000) * 0.065;
				total_bill += extra_charge;
			}
		}
	printf("\nThe total bill for account number %d is $%.2f\n\n",
		account_number, total_bill);
}



