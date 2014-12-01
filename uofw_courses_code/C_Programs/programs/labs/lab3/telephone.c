/* Author: Malachi Griffith
*  Date: Sept. 18 2002
*  Purpose: Computes telephone bills based on input from user
*/

#include <stdio.h>

main()
{
	double 	base_charge = 12.95,
		final_bill,
		bill,	
		tax,	
		long_charge,
		local_charge;
	int	local_calls;

	printf("Your base charge before taxes is %4.2f \n", base_charge);
	printf("Enter the number of local calls > ");
	scanf("%d", &local_calls); 
	printf("Enter the charge for long distance calls > ");
	scanf("%lf", &long_charge);

	local_charge = local_calls * 0.1;
	bill = base_charge + local_charge;
	bill += long_charge;
	tax = bill*.0872;
	bill *= 1.0872; 
	

	printf("The number of local calls was > %d\n", local_calls);	
	printf("The cost of local calls was > %5.2f\n", local_charge);
	printf("The cost of long distance calls was > %5.2f\n", long_charge);
	printf("The total tax was > %5.2f\n", tax);
	printf("The final bill is > %5.2f\n", bill);
} 




