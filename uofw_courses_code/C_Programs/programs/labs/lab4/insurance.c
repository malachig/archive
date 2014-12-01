/* Author: Malachi Griffith 
*  Date: Sept. 26 2002
*  Purpose: Illustrate the use of if statements to calculate 
*  insurance costs for different types of vehicles.
*/
#include <stdio.h>

main()
{
	int type;
	int class;
	double basic_fee = 30.0;
	double total_fee;

	printf("\n\nYour basic charge is $%4.2f\n\n", basic_fee);
	
	printf("Enter the type of vehicle, 1 for luxury, 2 for others > ");
	scanf(" %d", &type);
	
	printf("\nEnter the class of vehicle, 1 for city, 2 for rural > ");
	scanf(" %d", &class);

	if (type == 1)
	total_fee = 15.0 + 30.0;
		else
		total_fee = 30.0;

	if (class == 1)
	total_fee += 125.0;

	printf("\n\nThe total cost for you insurance is $%.2f\n\n", total_fee);
}
	
