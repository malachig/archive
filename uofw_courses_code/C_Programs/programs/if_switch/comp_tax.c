/* comp_tax.c */

/* Author: Malachi Griffith
*  Date: Oct. 2 2002
*  Purpose: Computes the tax due based on a tax table.
*  	    Pre: salary is defined
*	    Post: Returns the tax due for 0.0 <= salary <= 150,000.00;
*	    	  returns -1.0 if salary is outside the table range.
*/

#include <stdio.h>

main()
{
	double salary;
	double tax;
	printf("\n\nPlease enter your salary with 2 decimal points > ");
	scanf(" %lf", &salary);
	if (salary < 0.0)
		tax = -1.0;

	else if (salary < 15000.00)	/* First range */
		tax = 0.15 * salary;
	else if (salary < 30000.00)	/* Second range */
		tax = (salary - 15000.00) * 0.18 + 2250.00;
	else if (salary < 50000.00)	/* Third range */
		tax = (salary - 30000.00) * 0.22 + 5400.00;
	else if (salary < 80000.00)	/* Fourth range */
		tax = (salary - 50000.00) * 0.27 + 11000.00;
	else if (salary < 150000.00)	/* Fifth range */
		tax = (salary - 80000.00) * 0.33 + 21600.00;
	else
		tax = -1.0;
	printf("\n\nYour taxes are %.2f\n", tax);
	return(0);
}
