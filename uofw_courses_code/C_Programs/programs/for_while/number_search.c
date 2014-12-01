/* number_search.c */

/* Author: Malachi Griffith
*  Date: Oct. 25 2002
*  Purpose: Find smallest, largest and average of a set of N numbers.
*/

#include <stdio.h>
#include <math.h>

main()
{
	int N = 0;
	double number;
	double average;
	double largest;
	double smallest;
	double sum = 0;
	double range;
	double square, sum_square = 0;
	double std_dev;
	
	printf("\nEnter a number, 99 to quit > ");	
	scanf("%lf", &number);
	
	largest = number;
	smallest = number;
	N++;
		
	while(number != 99)
	{
		if (number > largest)
			largest = number;
		if (number < smallest)
			smallest = number;
		sum += number;
		square = pow(number, 2);
		sum_square += square;
	
		printf("Enter a number, 99 to quit > ");
		scanf("%lf", &number);
		N++;
	}
	average = sum / (N-1);
	std_dev = sqrt((sum_square / (N-1)) - pow(average, 2));
		
	printf("\nThe smallest number entered is: %.1f", smallest);	
	printf("\nThe largest number entered is: %.1f", largest);
	printf("\nThe average is: %.1f", average);
	printf("\nThe range of the values is: %.1f", largest - smallest);
	printf("\nThe standard deviation is: %.1f\n\n", std_dev);	
}



