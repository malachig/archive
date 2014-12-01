/* Author: Malachi Griffith
*  Date: Nov. 23 2002
*  Purpose: Reads characters from an input file and echos them to an
*  output file until EOF is reached.  At that time a summary of 
*  the number of upper case, lower case, exclamation points, and space
*  characters is printed.
*/

#include <stdio.h>

main()
{
	int spaces = 0, upper = 0, lower = 0, punc = 0, lines = 0;
	char ch;
	FILE *input_data;
	FILE *output_data;

	input_data = fopen("text.dat", "r");
	output_data = fopen("text.print", "w");

	fscanf(input_data, "%c", &ch);

	while (!feof(input_data))
	{
		if (islower(ch))
			lower++;
		else if (isupper(ch))
			upper++;
		else if (ispunct(ch))
			punc++;
		else if (isspace(ch))
			spaces++;
		
		if (ch == '\n')
			lines++;

	fprintf(output_data, "%c", ch);
	fscanf(input_data, "%c", &ch);
	}

	fprintf(output_data, "\n");
	fprintf(output_data, "Number of lines is %d\n", lines);
	fprintf(output_data, "Number of uppercase is %d\n", upper); 
	fprintf(output_data, "Number of lowercase is %d\n", lower);
	fprintf(output_data, "Number of punctuations is %d\n", punc); 

	fclose(input_data);
	fclose(output_data);
}



