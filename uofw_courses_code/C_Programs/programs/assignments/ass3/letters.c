/* Assignment 3: Question 1 (10 Marks) letters.c */

/* Author: Malachi Griffith
*  Date: Nov. 11 2002
*  Purpose: This program analyzes a data file containing two brief sentences.
*  It then prints the number of occurences of each lower and upper case 
*  letter of the alphabet.  Letters that do not appear in the data file will
*  not be shown as 0, but simply omitted from the print out.
*/

/* Note: Upper case letters correspond to 65-90 (A-Z) in ASCII.
*	 Lower case letters correspond to 97-122 (a-z) in ASCII.
*/


#include <stdio.h>
#define MAX_DATA 100   /* Max number of letters processed from data file */
#define MAX_LETTERS 26 /* Max array size for tallying number of each letter.*/

/* Function Prototype */
void read_array(char data[], int *size);
void process_array(const char data[], int size, int upper[], int lower[]);
void print_array(const int upper[], const int lower[]);

main()
{
	int size_local;
	char letters[MAX_DATA];
	
	int upper[MAX_LETTERS] = {0};
	int lower[MAX_LETTERS] = {0};

	read_array(letters, &size_local);

	process_array(letters, size_local, upper, lower);

	print_array(upper, lower);
 
	return(0);
}

/*
*  Function: read_array
*/
void 
read_array(char data[], int *size)
{
	char ch;
	int i = 0;

	/* Open the input file */
	FILE *input_data;
	input_data = fopen("letters.dat", "r");

	fscanf(input_data, "%c", &ch);

	while(!feof(input_data))
	{
		data[i] = ch;
		fscanf(input_data, "%c", &ch);
		i++;
	}
	*size = i;
	fclose(input_data);
}

/*
*  Function: process_array
*/
void 
process_array(const char data[], int size, int upper[], int lower[])
{
	int i;
	char ch_value;

	for(i = 0; i <= size; i++)
	{
		if (data[i] >= 65 && data[i] <= 90)
		{
			ch_value = data[i] - 65;
			upper[ch_value]++;
		}		
		
		if (data[i] >= 97 && data[i] <= 122)
		{
			ch_value = data[i] - 97;
			lower[ch_value]++;
		}
	}
}


/*
*  Function: print_array
*/
void 
print_array(const int upper[], const int lower[])
{
	char ch;
	int i;

	/* Open the output file */
	FILE *output_data;
	output_data = fopen("letters.print", "w");
	
	fprintf(output_data, "** The Occurence of Upper Case Letters **\n\n");

	ch = 'A';

	for(i = 0; i <= 25; i++)
	{
		if (upper[i] > 0)
			fprintf(output_data, 
				"\tThe occurence of letter %c was %d\n",
				ch, upper[i]);
		ch++;
	}
	
	fprintf(output_data,"\n** The Occurence of Lower Case Letters **\n\n");
	
	ch = 'a';

	for(i = 0; i <= 25; i++)
	{
		if (lower[i] > 0)
			fprintf(output_data,
				"\tThe occurence of letter %c was %d\n",
				ch, lower[i]);
		ch++;
	}

	fclose(output_data);
}


