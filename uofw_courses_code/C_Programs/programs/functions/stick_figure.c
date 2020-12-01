/* Author: Malachi Griffith
*  Date: Oct. 12 2002
*  Purpose: Uses a set of functions to draw a stick figure.
*/

#include <stdio.h>

/* Function prototypes */
void draw_circle(void);
void draw_intersect(void);
void draw_base(void);
void draw_triangle(void);
void draw_rectangle(void);
void skip_5_lines(void);
int
main(void)
{
	/* Draw a rocket */
	draw_triangle();
	draw_rectangle();
	draw_intersect();
	skip_5_lines();

	/* Draw a male stick figure */
	draw_circle();
	draw_rectangle();
	draw_intersect();
	skip_5_lines();

	/* Draw a female standing on a male */
	draw_circle();
	draw_triangle();
	draw_intersect();
	draw_circle();
	draw_rectangle();
	draw_intersect();
	skip_5_lines();

	printf("\n\n");
	return(0);
}

/*
*  Draws a circle 
*/

void
draw_circle(void)
{
	printf("    *   \n");
	printf("  *   * \n");
	printf("   * *  \n");
}

/*
*  Draws intersecting lines 
*/

void
draw_intersect(void)
{
	printf("   / \\  \n"); /* Use 2 \'s to print 1 */
	printf("  /   \\ \n");
	printf(" /     \\\n");
}

/*
* Draws a base line
*/

void
draw_base(void)
{
	printf(" -------\n");
}

/*
*  Draws a triangle
*/

void
draw_triangle(void)
{
	draw_intersect();
	draw_base();
}


/* Function rectangle */
void
draw_rectangle(void)
{
        printf(" -------\n");
        printf(" |     |\n");
        printf(" |     |\n");
        printf(" |     |\n");
        printf(" -------\n");
}
 
/* Function skip_5_lines */
void
skip_5_lines(void)
{
printf("\n\n\n\n\n");
}
