/* Author: Malachi Griffith
*  Date: Oct 12 2002
*  Purpose: Illustrates the use of simple functions
*/

#include <stdio.h>

/* Function prototypes */
void draw_rectangle(void);
void draw_parallel(void);
void draw_base(void);

main()
{
	printf("\n\n");
	draw_rectangle();
	printf("\n\n");
}

/*
* Draw Rectangle function calls on two other functions to 
* make the rectangle
*/

void
draw_rectangle(void)
{
	draw_base();
	draw_parallel();
	draw_base();
}

/*
* Draws a horizontal line
*/

void
draw_base(void)
{
printf("  ----------  \n");
}

/*
* Draws a pair of vertical lines
*/

void
draw_parallel(void)
{
printf(" |          |  \n");
printf(" |          |  \n");
printf(" |          |  \n");
printf(" |          |  \n");
printf(" |          |  \n");
printf(" |          |  \n");
printf(" |          |  \n");
printf(" |          |  \n");
printf(" |          |  \n");
printf(" |          |  \n");
}

