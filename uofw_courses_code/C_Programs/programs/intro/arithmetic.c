/*Arithmetic*/

/* Author: Malachi Griffith
*  Date: Sept. 16 2002 
*  Purpose: Simple math calculations */

#include <stdio.h>

main()
{
	/* Order of operations */

	printf("\n4*11/44 -19/2 is : %d\n", 4*11/44 - 19/2);
	printf("4*(11/44)-19/2 is : %d\n", 4*(11/44) - 19/2);


	/* Mixed versus integer division */
	
	printf("\ninteger division gives %d\n", 5/2);
	printf("mixed division gives %f\n", 5.0/2);


	/* Type cast - forces data to take another form */

	printf("\n%d\n", (int) 5.0/2);
	printf("%f\n", (float) 5/2);



}

