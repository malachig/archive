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
#define MIN_MARK 85

main()
{
	int mark[NUM_QUIZZES];
	int quiz;

	/* Read the data into the array "mark" */
	for (quiz = 0; quiz < NUM_QUIZZES; quiz++)
	{
		scanf("%d", &mark[quiz]);
		printf("Mark for Quiz %d = %d\n", quiz + 1, mark[quiz]);
	}

	for (quiz = 0; quiz < NUM_QUIZZES; quiz++)
		if (mark[quiz] >= MIN_MARK)
		{
			printf("The first mark of at least %d is: \n",
			       MIN_MARK);

			printf("Quiz #%d Mark: %3d\n", quiz + 1,
			       mark[quiz]);
			break;
		}
		
		if (quiz >= NUM_QUIZZES)
			printf("No mark greater than %d was found.\n",
			       MIN_MARK);
}
