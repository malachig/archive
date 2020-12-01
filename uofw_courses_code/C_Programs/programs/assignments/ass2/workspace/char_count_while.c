/* Assignment 2: Qustion 2 Part A. (10 Marks) char_count.c */

/* Author: Malachi Griffith
*  Date: Oct. 13 2002
*  Purpose: Use a batch program to read lines from a datafile called
*  	    "char_count.dat" and count the number of characters in each
*	    line.  Spaces and end-of-line marks will not be counted.
*/

#include <stdio.h>

main()
{
	char ch;
	int char_count = 0;
	int line_count = 0;

	FILE * input_data;

	input_data = fopen("char_count.dat", "r");

	fscanf(input_data, "%c", &ch);

	while (!feof(input_data))
	{
		while(ch != '\n')	
		{
			if (ch != ' ')
				char_count++;
			
			fscanf(input_data, "%c", &ch);	
		}
		
		line_count++;	
	
		if (char_count > 0)	
			printf("\nThe number of characters in line %d is %d",
		               line_count, char_count);
	
		fscanf(input_data, "%c", &ch);
		char_count = 0;
	} 	

	printf("\n\n");	
	fclose(input_data);	
	return(0);
}



