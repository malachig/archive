/* Author: Malachi Griffith
*  Date: Nov. 17 2002
*  Purpose: Illustrate sting input/output with printf and scanf
*/

#include <stdio.h>

#define STRING_LEN 10

int
main(void)
{
	char dept[STRING_LEN];
	int course_num;
	char days[STRING_LEN];
	int time;

	printf("Enter department code, course number, days and ");
	printf("time like this: \n> COSC 2060 MWF 1410\n> ");
	scanf("%s%d%s%d", dept, &course_num, days, &time);
	printf("%s %d meets %s at %d\n", dept, course_num, days, time);

	return(0);	
}



