/* Author: Malachi Griffith
*  Date: Nov. 2 2002
*  Purpose: Illustrate the basic use of pointers.
*/

#include <stdio.h>

void test (float *);

main()
{
	float mark = 78.3;

	test(&mark);

	printf("In main, mark is %.2f\n", mark);
}


void test (float *mark_ptr)
{
	*mark_ptr = 40.5;
	printf("In test, mark is %.2f\n", *mark_ptr);
}


