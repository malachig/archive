/*degree_convert.c*/
/* Author: Malachi Griffith
*  Date:  Sept. 20 2002 
*  Purpose:
*/

#include <stdio.h>
main()
{
	double farenheit,
	       celsius;
	printf("\nType in the temperature in degrees farenheit >");
	scanf("%lf", &farenheit);
	celsius = ((5.0/9.0) * (farenheit - 32));
	printf("\nThe corresponding temperature in degrees celsius is %6.2f\n\n",
	celsius);

}
