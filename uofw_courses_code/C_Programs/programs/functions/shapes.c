/*
*  Author: Malachi Griffith
*  Date: Oct. 21 2002
*  Purpose: One function draws a triangle the other draws a rectangle
*/

#include <stdio.h>

/* Function prototypes */
void rectangle(void);
void triangle(void);

main()
{
	triangle();
	rectangle();
	printf("\n\n");
	triangle();
	rectangle();

	printf("\n\n");
	return(0);

}

/* Function rectangle */
void
rectangle(void)
{
	printf("\n\t*****");
	printf("\n\t*   *");
	printf("\n\t*   *");
	printf("\n\t*****");
}

/* Function triangle */
void
triangle(void)
{
	printf("\n\t  *  ");
	printf("\n\t * * ");
	printf("\n\t*****");
}


