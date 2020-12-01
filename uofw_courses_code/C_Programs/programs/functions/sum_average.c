/* Author: Malachi Griffith
*  Date: Oct. 8 2002 
*  Purpose: Compute the sum and average of two numbers. 
*/

#include <stdio.h>

int
main(void)
{
	int    one, 		/* input number 1 */
	       two,		/* input number 2 */
	       sum,		/* Compute average of the numbers */
	       average; 	/* Display sum and average. */

	/* Get numbers from user */
	printf("\nPlease enter two integer values seperated by a space > ");
	scanf("%d%d", &one, &two);

	/* Calculate average and sum */
	sum = one + two;
	average = sum / 2;

	printf("The sum and average of these numbers is: %d and %d\n:", sum, average);
	
	return(0);
}



