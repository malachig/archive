/* Author: Malachi Griffith
*  Date: Nov. 3 2002  
*  Purpose: Simple function using return statment.
*/

#include <stdio.h>

/* function prototype */
int sum(int a, int b);

int 
main(void)
{
	int x, y, z;
	
	x = 5; y = 3;

	printf("   x   y   z\n\n");

	z = sum(x, y);
	printf("%4d%4d%4d\n", x, y, z);

	z = sum(y, x);
	printf("%4d%4d%4d\n", x, y, z);

	x = sum(z, y); 
	printf("%4d%4d%4d\n", x, y, z);

	x = sum(z, z);
	printf("%4d%4d%4d\n", x, y, z);
	
	y = sum(y, y);
	printf("%4d%4d%4d\n", x, y, z);

	return(0);
}

int
sum(int a, int b)
{
	int temp;
	temp = a + b;
	return(temp);
}


