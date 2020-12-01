/* Author: Malachi Griffith
*  Date: Sept. 18 2002
*  Purpose:  Simple calculation of salary, correcting for a yearly increase
*/

#include <stdio.h>

main()
{
	double salary,
	       increase,
	       year1,
	       year2,
	       year3;
		
	printf("\nEnter your annual salary as of the current year > ");
	scanf("%lf", &salary);
	printf("Enter the estimated yearly rate of growth in %");
	scanf("%lf", &increase);
	
	year1 = salary + (salary * increase/100.0);
	year2 = year1 + (year1 * increase/100.0);
	year3 = year2 + (year2 * increase/100.0);

	printf("Your salary for next year is > %7.2f \n", year1);
	printf("Your salary for the following year is > %7.2f \n", year2); 
	printf("Your salary for the third year is > %7.2f \n\n", year3);
}
  	



