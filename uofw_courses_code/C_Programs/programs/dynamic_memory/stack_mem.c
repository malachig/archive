/* Author: Malachi Griffith
*  Date: Nov. 21 2002
*  Purpose: Illustrate the use of Stack memory.  Stack memory is used 
*  to store variables in functions because they may be called multiple
*  times.
*/

#include <stdio.h>
#define SIZE_CONSTANT 10

int iglobal = 2;

/* Function Prototype */
int func(int []);

main()
{
	int imain[SIZE_CONSTANT] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	int ireturn;

	ireturn = func(imain);
	printf("ireturn = %d in main\n", ireturn);
}

int
func(int ipass[SIZE_CONSTANT])
{
	int ifunc;
	int i;
	int sum = 0;
		
	for(i = 0; i < SIZE_CONSTANT; i++)
		sum += ipass[i];
	
	ifunc = sum/iglobal;
	return (ifunc);
}



