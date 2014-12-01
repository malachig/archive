/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: A program for reading in quiz marks into a one-dimensional 
*  array.  When the input cursor appears you have to enter in 10 data values
*  then press <ctrl> d and <enter>, to end the input and continue with 
*  the program.
*/

#include <stdio.h>

#define NUM_QUIZZES 10

main()
{
	int mark[NUM_QUIZZES];
	int quiz;
	int sum = 0;

	double average;

	for (quiz = 0; quiz < NUM_QUIZZES; quiz++)
	{
		scanf("%d", &mark[quiz]);
		printf("Mark for Quiz %d = %d\n", quiz + 1, mark[quiz]);
		sum += mark[quiz];
	}

	average = (double) sum / NUM_QUIZZES;

	printf("The average quiz mark is: %.1f\n", average);
}



