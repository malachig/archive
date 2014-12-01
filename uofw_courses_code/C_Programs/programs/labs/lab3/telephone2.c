/* Author: Malachi Griffith
*  Date: Sept. 18 2002
*  Purpose: Computes telephone bills based on input from user
*  The program also illustrates the use of assignment operators such as "+="
*  These operators are used as the "=" sign might be but instead of simply 
*  assigning a value directly they increase, decrease, multiply or divide the value.
*/

#include <stdio.h>

main()
{
	double 	base_charge = 12.95,
		bill = 0,	
		tax,	
		long_charge,
		local_charge;
	int	local_calls;

	printf("\nYour base charge before taxes is %4.2f \n", base_charge);
	printf("Enter the number of local calls > ");
	scanf("%d", &local_calls); 
	printf("Enter the charge for long distance calls > \n");
	scanf("%lf", &long_charge);

	printf("the number of local calls was > %d\n", local_calls);
	local_charge = local_calls * 0.1;
	printf("The cost of the local calls was > %5.2f\n", local_charge);
	printf("The cost of long distance calls was > %5.2f\n", long_charge);
	
	bill += base_charge;
	bill += long_charge;
	bill += local_charge;
	printf("The total bill before taxes is > %5.2f\n", bill);
	tax = bill*.0872;
	printf("The taxes on this bill are > %5.2f\n", tax);
	bill *= 1.0872; 
	printf("Your final bill including taxes is > %5.2f\n", bill);

} 
