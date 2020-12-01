/* Author: Malachi Griffith
*  Date: Nov. 2 2002 
*  Purpose: Illustrate the use of simple arguments in functions (not pointers).
*/

#include <stdio.h>

void test(float);

main()
{
	float mark = 78.3;
	test(67.4);
	printf("In main, mark is %.2f\n", mark);
}

void
test(float mark)
{
	printf("In test, mark is %.2f\n", mark);
}
