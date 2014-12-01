/* Author: Malachi Griffith
*  Date: Nov. 6 2002
*  Purpose: Illustrate how hexadecimel counting works.  The hexadecimel 
*  system is made up of 16 digits (0,1,2,3,4,5,6,7,8,9,a,b,c,d,e,f).
*  The address of a pointer will be printed and then incremented by one. 
*  For an integer this means that the adress is incremented by four bytes,
*  the amount of space required to store one integer.
*/

#include <stdio.h>

main()
{
	int *value;
	int i;
	
	for(i = 0; i<=30; i++)
	{
		printf("\nThe address is > %p", value);
		value += 1;
	}

	printf("\n");
}



