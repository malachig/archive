/* Author: Malachi Griffith
*  Date: Sept.  20 2002
*  Purpose: Illustrate the basic usage and syntax of the if\else selection structure.
*/

#include <stdio.h>

main()
{

	double grade;
	printf("\n\nEnter the grade you recieved in the course > ");
	scanf(" %lf", &grade);
	if (grade >= 50)
		printf ("passed\n\n");
	else
		printf("Failed\n\n");

}
