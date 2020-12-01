/* while.c */

/* Author: Malachi Griffith
*  Date: 
*  Purpose:
*/

#include <stdio.h>

main()
{
	int exams_marked,
	    students;

	exams_marked = 0;
	printf("\n\nEnter the number of students who wrote the exam > ");
	scanf(" %d", &students);
	
	while (exams_marked <= students)
		{
		printf("Exams marked = %d\n", exams_marked);
		exams_marked++;
		}
}



