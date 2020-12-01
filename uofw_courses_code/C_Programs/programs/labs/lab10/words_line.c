/* words_line.c */

/* Author: Malachi Griffith
*  Date: Nov. 7 2002 
*  Purpose: Will count the number of words in a line.
*/

#include <stdio.h>

#define MAX_SIZE 81

main()
{
	char string[MAX_SIZE];
	int i;
	int word_count = 0;

	FILE *input_string;
	input_string = fopen("sentence.dat", "r");

	fgets(string, MAX_SIZE, input_string);
	
	printf("\n");

	for(i = 0; string[i] != '\0'; i++)
	{
		printf("%c", string[i]);
		if (string[i] == ' ' || string[i] == '.')
			word_count++;	
	}
	printf("\nThe word count for this sentence is %d\n\n", word_count);
}





