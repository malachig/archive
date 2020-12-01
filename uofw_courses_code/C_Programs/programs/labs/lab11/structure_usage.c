/* structure_usage.c */

/* Author: Malachi Griffith
*  Date: Nov. 14 2002 
*  Purpose: Illsutrate the usage of some simple structures
*/

#include <stdio.h>

/* Define a structure type which contains: a 30-character name, an integer
 * student ID number, a 5-char year ID, and an array of five 15-character 
 * course names.
*/

/* First by using a structure tagname, Student_tag. */

	struct 	Student_tag{
		char student_name[30];
		int student_id;
		char year_id[5];
		char courses[5][15];
		}
/* Second by using a typedef statement */

	typedef struct{
		char student_name[30];
		int student_id;
		char year_id;
		char courses[5][15];
		} Student_Struct;

main()
{
	struct Student_tag student1;

	Student_Struct student2;

}



