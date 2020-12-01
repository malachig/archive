/*number_format.c*/
/* Author: Malachi Griffith
*  Date: Sept. 18 2002
*  Purpose: To illustrate how to format numbers for program output
*/

#include <stdio.h>
 
main()
{
int a = 504;
float b = 302.558;
float c = -12.31;

printf("%5d", a);
printf("%11.2f", b);
printf("%9.1f\n", c);
}

