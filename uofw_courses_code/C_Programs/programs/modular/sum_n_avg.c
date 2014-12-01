/* Author: Malachi Griffith
*  Date: Nov. 3 2002  
*  Purpose: Use a function with pointers to calculate the sum and average 
*  of three variables and return them to main.
*/

#include <stdio.h>

/* Function Prototypes */
void sum_n_avg (double x, double y, double z, double *sum, double *avg);

main()
{
	double x, y, z;
	double sum1, average;

	printf("\nPlease enter three numbers seperated by spaces > ");
	scanf("%lf%lf%lf", &x, &y, &z);

	sum_n_avg(x, y, z, &sum1, &average); 

	printf("\nThe sum is %.1f and the average is %.1f", sum1, average);
	printf("\n\n");

}

void
sum_n_avg (double x, double y, double z, double *sum, double *avg)
{

	*sum = x + y + z;
	*avg = *sum / 3;
}

