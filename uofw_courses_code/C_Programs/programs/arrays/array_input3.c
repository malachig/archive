/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: A program for reading in quiz marks into a one-dimensional 
*  array.  When the input cursor appears you have to enter in 10 data values
*  then press <ctrl> d and <enter>, to end the input and continue with 
*  the program.
*  Also finds the highest data value, the element containing this mark,
*  and the average of all marks.
*/

#include <stdio.h>

#define NUM_QUIZZES 10

main()
{
	int mark[NUM_QUIZZES];

	int quiz;
	int best_quiz;
	int sum = 0;
	int max_mark;

	double average;

	/* Read the data into the array "mark" */
	for (quiz = 0; quiz < NUM_QUIZZES; quiz++)
	{
		scanf("%d", &mark[quiz]);
		printf("Mark for Quiz %d = %d\n", quiz + 1, mark[quiz]);
		sum += mark[quiz];
	}

	max_mark = mark[0];
	best_quiz = 0;

	for (quiz = 0; quiz < NUM_QUIZZES; quiz++)
	{
		sum += mark[quiz];
	
		if (mark[quiz] > max_mark)
		{
			max_mark = mark[quiz];
			best_quiz = quiz;
		}
	} 

	average = (double) sum / NUM_QUIZZES;

	printf("The average quiz mark is: %.1f\n", average);
	printf("The highest mark was on quiz %d\n", best_quiz+1);
	printf("The highest mark was %d\n", max_mark);
}
