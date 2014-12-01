/* coin_toss.c */

/* Author: Malachi Griffith
*  Date: Nov. 13 2002
*  Purpose: Read lines from a data file, then use a random number
*/

#include <stdio.h>
#include <string.h> /* needed for string functions */
#include <stdlib.h> /* needed for srand() */
#include <time.h>   /* needed for time() */

#define MAX_LENGTH 81
	
	typedef enum {heads, tails} outcome_t;

main()
{
	outcome_t coin_toss;
	int random_num;
	
	char input_str[MAX_LENGTH];
	
	/* Open input and output files */
	FILE *input_data;
	FILE *output_data2;
	FILE *output_data1;
	
	input_data = fopen("coin_toss.dat", "r");
	output_data1 = fopen("a_file.dat", "w");
	output_data2 = fopen("b_file.dat", "w");
	
	fprintf(output_data1, "\n");
	fprintf(output_data2, "\n");
			
	fgets(input_str, MAX_LENGTH, input_data);

	srand(time(0));

	while (!feof(input_data))
	{
		random_num = rand();

		if ((random_num % 2) == 0)
			coin_toss = heads;
		else 
			coin_toss = tails;

		if (coin_toss == heads)
			fputs(input_str, output_data1);
		
		else if (coin_toss == tails)
			fputs(input_str, output_data2);
	
		fgets(input_str, MAX_LENGTH, input_data);
	}
	
	fprintf(output_data1, "\n");
	fprintf(output_data2, "\n");

	/* Close the data files */	
	fclose(input_data);
	fclose(output_data1);
	fclose(output_data2);
}



