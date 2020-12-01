/*  Assignment 3: Question 2 (15 Marks) sort_names.c */

/* Author: Malachi Griffith
*  Date: Nov. 11 2002 
*  Purpose: This program will sort the names provided in the data file
*  "assign3.dat".  Each line contains the surname followed by the first
*  name of a graduating member of the class.  Each line is at most 25 
*  characters long, and there at most 20 students.  The actual size is 
*  known only after all data is read.
*/

#include <stdio.h>
#include <string.h>

#define MAX_LETTERS 25
#define MAX_STUDENTS 20

/* Function Prototypes */
int read_array(char names[][]);
void sort_array(char names[][], int num_students);
void display_array(char names[][], int num_students);
int find_min(int start, int end, char names[][]);

main()
{
	int students;  /* Actual number of students found in datafile */ 

	char student_names[MAX_STUDENTS][MAX_LETTERS];

	students = read_array(student_names);

	sort_array(student_names, students); 

	display_array(student_names, students); 

	return(0);
}

/*
*  Function: read_array
*/
int 
read_array(char names[][MAX_LETTERS])
{
	int current = 0; /* Control variable for loop, increases after
			  * each student name is processed. */

	/* Open the input file */	
	FILE *input_data;
	input_data = fopen("assign3.dat", "r");
	
	while(!feof(input_data))
	{
		fgets(names[current], MAX_LETTERS, input_data);	

		current++;
	}
	return(current - 1);	
}

/*
*  Function: sort_array - uses selection sort and helper function (find_min)
*/
void 
sort_array(char names[][MAX_LETTERS], int num_students)
{
	int current_min;
	int index;
	char temp_array[MAX_LETTERS];

	for(index = 0; index < num_students; index++)
	{
		current_min = find_min(index, num_students, 
					names); 
			
		if(index != current_min)
		{

			strcpy(temp_array, names[index]);
	
			strcpy(names[index], names[current_min]);
		
			strcpy(names[current_min], temp_array);
		
		}
	}
}

/*
*  Function: find_min
*/
int 
find_min(int start, int end, char names[][MAX_LETTERS])
{
	int current_min;
	int index;
	int compare_string;

	current_min = start;
	
	for(index = start + 1; index < end; index++)
		{
		/* Use string compare to compare each string to the starting
		*  string (current).  If the test string is smaller than
		*  the string tested a -ve result is returned and this means
		*  that the test string comes first alphabetically.  The 
		*  position of this string will be considered the minimum
		*  until an even smaller one is found.*/
	
 
		compare_string = strcmp(names[index], names[current_min]);
		if (compare_string < 0)
			current_min = index;
		}
	return(current_min);
}
	
void 
display_array(char names[][MAX_LETTERS], int num_students)
{
	int student;

	/* Open file for print out */	
	FILE *output_data;	
	output_data = fopen("sort_names.print", "w");

	for (student = 0; student < num_students; student++)
		fputs(names[student], output_data);	
	
	fclose(output_data);
}
