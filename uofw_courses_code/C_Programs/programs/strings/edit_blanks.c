/* Author: Malachi Griffith
*  Date: Nov. 17 2002 
*  Purpose: The function trim blanks takes as it's input a string with leading
*  and trailing blanks.  It then removes the leadingand trailing blanks and 
*  returns the edited string.
*/

#include <stdio.h>
#include <string.h>

#define MAX_LENGTH 81

/* Function Prototype */
void trim_blanks(const char to_trim[], char trimmed[]);

main()
{
	char to_trim[MAX_LENGTH] = {"Hello this is me"};
	char trimmed[MAX_LENGTH];
	
	printf("Enter a string please, \n>");
	gets(to_trim);

	trim_blanks(to_trim, trimmed);

	printf("\nThe edited string is: \n");
	puts(trimmed);
}

/*
*  Function: trim_blanks
*/
void
trim_blanks(const char to_trim[], char trimmed[])
{
	int i;
	int first_non_blank;
	int last_non_blank;
	int difference;
	
	/* Find subscript of first nonblank in to_trim */
	for(i = 0; i <= MAX_LENGTH - 1; i++)
	{
		if(to_trim[i] != ' ')
		{
			first_non_blank = i;
			break;
		}
	}

	/* find subscript of last nonblank in to_trim */
	for(i = first_non_blank; i <= MAX_LENGTH - 1; i++)
	{
		if(to_trim[i] != ' ' && to_trim[i] != '\n' && to_trim[i] != '\0')
			last_non_blank = i;
	}
	difference = (last_non_blank - first_non_blank) + 1;

	/* Use strncpy to store trimmed string in trimmed */
		strncpy(trimmed, &to_trim[first_non_blank], difference);
		trimmed[difference] = '\0';
}



