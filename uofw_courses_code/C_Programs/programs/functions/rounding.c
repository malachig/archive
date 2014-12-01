/* Author: Malachi Griffith
*  Date: Oct. 22 2002  
*  Purpose: Provide a function that rounds numbers to 2 decimel points
*  Use the math function 'floor(x)'.  This function returns the largest
*  integral value that is not larger than x.
*/

#include <stdio.h>
#include <math.h>

/* Function prototype */
double round(int num_dec, double number);

main()
{
	double final;
	int decimels;
	double x;

	printf("\nEnter the value you wish to round > ");
	scanf("%lf", &x);
	printf("How many decimel points would you like it rounded to? ");
	scanf("%d", &decimels);

	final = round(decimels, x);

	printf("\n\nThe rounded value is: %f\n\n", final);
}



/*
*  Function: round
*  Pre: the number of decimel places to round to as well as the actual
*  number to be rounded are defined.
*/

double
round(int num_dec, double number) 
{

	double trunc_x;
	double x_pow;
	double trunc_x_pow;
	double final_x;
	int power;

	power = pow(10, num_dec);

	trunc_x = number - floor(number);
	x_pow = trunc_x * power;
		
	trunc_x_pow = x_pow - floor(x_pow);

	if (trunc_x_pow >= 0.5)
	{
		final_x = floor(x_pow) + 1;
		final_x /= power;
		final_x += floor(number);
	}

	else
	{
		final_x = floor(x_pow);
		final_x /= power;
		final_x += floor(number);
	}

	return(final_x);
}



