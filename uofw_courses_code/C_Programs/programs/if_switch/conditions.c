/* conditions.c */

/* Author: Malachi Griffith
*  Date: Sept. 21 2002
*  Purpose: To illustrate the use of simple conditional statements.
*/

#include <stdio.h>

main()
{
	int year, even;
	double test;
	double age, 
	       water,
	       speed,
	       y = 3.0,
	       w,
	       x = 2.0,
	       z = 5.0;
 	printf("\n\nEnter a value for age > ");
	scanf(" %lf", &age);
	test = (age >= 18 && 21 >= age);
	printf("Test of whether age is between 18 and 21: %1.0f\n\n", test);

	printf("Enter a value for water > ");
	scanf(" %lf", &water);
	test = (water < 1.5 && 0.1 < water);
	printf("Test of whether water is between 0.1 and 1.5: %1.0f\n\n", test);

	printf("Enter the year > ");
	scanf(" %d", &year);
	even = (year % 4 == 0);
	printf("Test of whether the year is divisible by 4: %d\n\n", even);

	
	printf("Enter the speed > ");
	scanf(" %lf", &speed);
	test = !(speed > 55);
	printf("Test of whether speed is not greater than 55: %1.0lf\n\n", test);

	test = (y >= x && z >= y);
	printf("Test of whether y is greater than x and less than z, x=2, y=3, z=5: %1.0f\n\n",
	test);

	printf("Enter a value of w > ");
	scanf(" %lf", &w);	
	test = (w == 6 || !(w > 3));
	printf("Test of whether w is either equal to 6 or not greater than 3: %1.0f\n\n", test);
}
