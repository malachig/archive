/* Author: Malachi Griffith
*  Date: Oct. 6 2002
*  Purpose: Compute the sum of a list of exam scores stored in 
*	     the file scores.dat. Using the EOF method of terminating a
*	     loop.   
*/

#include <stdio.h>

int
main(void)
{
	FILE * ifptr;	/* Input file pointer */
	
	
	int sum = 0,	/* Output - sum of scores input so far */
	    score,	/* Input - current score */
	    input_status; /* Status value returned by fscanf */
	ifptr = fopen("scores.dat", "r");

	printf("Scores\n");

	input_status = fscanf(ifptr, "%d", &score);
	while (input_status != EOF)	
	/* Could also use 'while(!feof(ifptr))' */	
	{
		printf("%5d\n", score);	
		sum += score;
		input_status = fscanf(ifptr, "%d", &score);
	}

	printf("\nSum of exam scores is %d\n", sum);
	fclose(ifptr);

	return(0);
}



