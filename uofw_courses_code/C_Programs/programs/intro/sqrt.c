/*sqrt.c*/
/* Author: Malachi Griffith 
* Date: Sept. 16 2002
* Purpose: Simply illustrates use of a square root function
* Note: Must be compiled with the statement "cc -lm sqrt.c"
* because it calls on the math.h library */

#include <stdio.h>
#include <math.h>


main()
{
double number1, a;
a = 9;
number1 = sqrt(a);
printf("The squareroot of 9 is %f. \n", number1);

return(0);
}
 
