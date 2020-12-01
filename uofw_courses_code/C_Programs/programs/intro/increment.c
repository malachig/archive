/*increment.c*/
/* Author: Malachi Griffith
*  Date: Sept. 16 2002
*  Purpose: Increment and Decrement functions
*/
main()
{
	int i = 5;
	float x = 2.5;

	printf("Original values: i = %d, x = %.2f\n", i, x);
	printf("In postfix position, i = %d, x = %.2f\n", i++, x--);
	
	printf("In 2nd printf after postfix, i=%d, x=%.2f\n", i, x);
	
	i = 5;
	x = 2.5;

	printf("In prefix position, i = %d, x = %.2f\n", ++i, --x);
 

}

