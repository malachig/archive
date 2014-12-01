/*  Assignment 3: Question 3 (15 Marks) edit_names.c */

/* Author: Malachi Griffith
*  Date: Nov. 11 2002 
*  Purpose: This program will edit the names provided in the data file
*  "sort_names.print".  Each line contains the surname followed by the first
*  name of a graduating member of the class.  Each line is at most 25 
*  characters long, and there at most 20 students.  The actual size is 
*  known only after all data is read.  This program will edit the names 
*  to fit the format "Smith F." and display them in an output file named
*  "edit_names.print".
*/

#include <stdio.h>
#include <string.h>

#define MAX_LETTERS 25
#define MAX_STUDENTS 20

/* Function Prototypes */
int find_edit_length(int length[]);
void edit_names(int num_students, const int length[], char names[][]);
void display_array(char names[][], int num_students);

main()
{
	int students;  /* Actual number of students found in datafile */ 

	int edit_length[MAX_LETTERS]; /* Number of characters up to and 
				       * including the first letter of the
				       * given name for each student */

	char student_names[MAX_STUDENTS][MAX_LETTERS];
	
	/* Function Calls */
	students = find_edit_length(edit_length);
	edit_names(students, edit_length, student_names);
	display_array(student_names, students); 

	return(0);
}

/*
*  Function: find_edit_length
*/
int
find_edit_length(int length[])
{
	int current_name = 0;
	int counter = 0;
	char ch;
	
	FILE *input_data;
	input_data = fopen("sort_names.print", "r");

	fscanf(input_data, "%c", &ch);
	
	while (!feof(input_data))
	{
		while (ch != ' ')
		{
			counter++;
			fscanf(input_data, "%c", &ch);
		}
	
		/* Now increase the counter to include the space and 
		 * next character */
		counter += 2;
	
		/* Place this value of counter in the length array.
		 * This will be the number of letter to process for the
		 * current name being read. */
		length[current_name] = counter;

		/* Now advance scanf cursor to end of line */
		while (ch != '\n')
			fscanf(input_data, "%c", &ch);
	
		current_name++;
		counter = 0;
		
		fscanf(input_data, "%c", &ch);
	}
	return(current_name);
}
		
/*
*  Function: read_array
*/
void
edit_names(int num_students, const int length[], char names[][MAX_LETTERS])
{
	int student;
	char junk_array[MAX_STUDENTS][MAX_LETTERS];

	/* Open the input file */	
	FILE *input_data;
	input_data = fopen("sort_names.print", "r");
	
	for(student = 0; student < num_students; student++)
	{
		fgets(names[student], length[student] + 1, input_data);	
		
		strcat(names[student], ".\n");

		fgets(junk_array[student], MAX_LETTERS, input_data);
	}
}

/*
*  Function: display_array
*/
void 
display_array(char names[][MAX_LETTERS], int num_students)
{
	int student;

	/* Open file for print out */	
	FILE *output_data;	
	output_data = fopen("edit_names.print", "w");

	for (student = 0; student < num_students; student++)
		fputs(names[student], output_data);	
	
	fclose(output_data);
}
