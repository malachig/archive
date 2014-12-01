/*grades.c*/
/* Author: Malachi Griffith
*  Date: Sept. 20 2002 
*  Purpose:
*/

#include <stdio.h>

main()
{
	char desired_grade;
	double	minimum_average,
		current_average,
		percent_counts,
		marks_needed,
		score_needed;
	
	printf("\n\nEnter the grade you desire > ");
	scanf(" %c", &desired_grade);

	printf("Enter the minimum mark required for this grade > ");
	scanf(" %lf", &minimum_average);

	printf("Enter your current average > ");
	scanf(" %lf", &current_average);
	 
	printf("Enter how much the final counts as a percentage of your grade > ");
	scanf(" %lf", &percent_counts);

	marks_needed = (double) (minimum_average - current_average);
	score_needed = (double) (marks_needed/(percent_counts/100.0));

	printf("\nThe score you need on the final test is: %4.2f\n\n", score_needed);

}
