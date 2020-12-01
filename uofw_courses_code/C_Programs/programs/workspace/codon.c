/*
* Author: Malachi Griffith
* Date: Sept. 14 2002
* Purpose: Illustrate use of scanf to take a series of letters
*         as input data and assign each to a seperate char variable.
*/

#include <stdio.h>

main()
{
	char first,
     	second,
     	third;
	char answer;
	char repeat = 1;

while (repeat == 1)
	{

	printf("\nEnter a three letter codon> ");
	scanf(" %c %c %c",&first,&second,&third);

	if (first == 'a' && second == 'a' && third == 'a')
		printf("\n  The corresponding amino acid is> Lysine\n\n");
	else if (first == 'a' && second == 't' && third == 'g')
		printf("\n  The corresponding amino acid is> Methionine\n\n");
	
	printf("Would you like to determine another amino acid?");
	scanf(" %c", &answer);
	if (answer == 'y') repeat = 1;
	else repeat = 0;
	}
	
	return(0);
}





