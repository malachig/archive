/* Author: Malachi Griffith
*  Date: Oct. 5 2002
*  Purpose: To display a table of degree values with their corresponding
*            sin and cosine values.
*/

#include <stdio.h>
#include <math.h>
main()
{
	int init_degree; 	/* Start of range chosen by user */
	int final_degree;	/* End of range chosen by user */
	int step_degree;	/* Increment for table */

	int value;		/* Control variable for for statement */

	double sin_value;
	double cos_value;

	printf("\n\nPlease enter the start of the range in degrees > ");
	scanf("%d", &init_degree);
	printf("Please enter the end of the range in degrees > ");
	scanf("%d", &final_degree);
	printf("Please enter the increment in degrees > ");
	scanf("%d", &step_degree);

	printf("\n\n\tDegree\t\tSin Value\t\tCos Value");

	for (value = init_degree;
	     value <= final_degree;
	     value += step_degree)
	{
	sin_value = sin ((double) (value / 57.295));
	cos_value = cos ((double) (value / 57.295));
	printf("\n\t%d\t\t%.4f\t\t%.4f", value, sin_value, cos_value);
	}
	printf("\n\n");
}



