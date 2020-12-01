/* Author: Malachi Griffith
*  Date: Nov. 2 2002
*  Purpose: Illustrate the basic use of pointers.
*  Call by reference allows us to affect variables in the calling function
*  in other ways than the return statement.
*  This program does not use pointers but the next one will!
*/

#include <stdio.h>

void swap (int, int); 

main()
{
	int a = 5;
	int b = 9;
	printf("Before swapping, in main, a = %d, b = %d\n", a, b);

	swap (a, b);

	printf("After swapping, in main, a = %d, b = %d\n", a, b);
}

void
swap(int a, int b)
{
	int temp;
	
	temp = a;
	a = b;
	b = temp;

	printf("In swap, a = %d, b = %d\n", a, b);
}


