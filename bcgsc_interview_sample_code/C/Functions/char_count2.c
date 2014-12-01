/* Assignment 2: Qustion 2 Part B. (5 Marks) char_count2.c */

/* Author: Malachi Griffith
*  Date: Oct. 13 2002
*  Purpose: Use a batch program to read lines from a datafile called
*  	    "char_count.dat" and count the number of characters in each
*	    line.  Spaces and end-of-line marks will not be counted.  
*	    The actual counting in this program is done by a sub-function
*	    called "countchar".
*/

#include <stdio.h>

/* Function prototype */
int countchar(int char_count, char ch);

main()
{
	char ch;		/* Stores characters from the datafile */
	int char_count = 0;	/* Counts number of characters on each line */
	int line_count = 0;     /* Keeps track of the line being counted */

	/* Define the datafiles */ 
	FILE *input_data;
	FILE *output_data;

	/* Open a file and get ready to read it */
	input_data = fopen("char_count.dat", "r");

	/* Open a file and get ready to write to it */
	output_data = fopen("char_count2.print", "w");

	/* Use the feof command to stop reading when the file ends */
	while (!feof(input_data))
	{
		/* The following loop continues to read characters and 
		*  count them (if they are not a space character) until
		*  an end-of-line character is reached.  Then the loop
		*  breaks and the rest of the while loop is executed. 
		*/
	
		for (fscanf(input_data, "%c", &ch);
		     (ch != '\n');
		     fscanf(input_data, "%c", &ch))	
		{
			/* Call the user defined function - char_count */
			char_count = countchar(char_count, ch);	
		}
		
		line_count++;	

		/* In case a line has no characters it will still be counted
		   but there is no point in displaying the value of 0.
		*/	
		if (char_count > 0)	
		{
			fprintf(output_data, "\nThe number of characters ");
		        fprintf(output_data, "in line %d is %d", line_count,
			       char_count);
		}
	
		char_count = 0; /* This variable must be reinitialised or
				*  the count accumulates from line to line
				   as the loop is repeated */
	} 	

	fprintf(output_data, "\n\n");	

	/* Close the data files after completion */
	fclose(input_data);
	fclose(output_data);

	return(0);
}

/* Function countchar */
/* This function simply receives characters from the function main and
*  if they are not a space, increases the value of the "char_count"
*  variable.
*  Pre: char_count and ch are defined.
*  Post: The updated value of char_count is returned to the main function. 
*/ 

int
countchar(int char_count, char ch)
{
	if (ch != ' ')
		char_count++;
	return(char_count);
}
