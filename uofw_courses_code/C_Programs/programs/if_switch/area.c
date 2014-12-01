/* area.c */

/* Author: Malachi Griffith
*  Date:  Sept. 28 2002
*  Purpose:  To calculate the area of a square or triangle
*/

#include <stdio.h>

main()
{
	char choice;
	double 	area = 0.0,
		base,
		height,
		side;

	printf("\n\nWould you like to find the area of a square or triangle, s or t > ");
	scanf(" %c", &choice);

	if (choice == 's' || choice == 'S')
		{
		printf("\nEnter the length of a side for the square > ");
		scanf(" %lf", &side);
		area = side * side;
		printf("\nThe area of this square is: %.2f\n\n", area);
		}
	if (choice == 't' || choice == 'T')	
		{
		printf("\nEnter the base of the triangle > ");
		scanf(" %lf", &base);
		printf("Enter the height of the triangle > ");
		scanf(" %lf", &height);

		area = (base * height) / 2;
		printf("\nThe area of this triangle is: %.2f\n\n", area);
		} 
		
	if (!(choice == 's' || choice == 'S' || choice == 't' || choice == 'T'))	
		printf("\nThat is not a valid entry!\n\n");
}



