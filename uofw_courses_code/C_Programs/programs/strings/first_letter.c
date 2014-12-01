/* Author: Malachi Griffith
*  Date: Nov. 17 2002 
*  Purpose: Takes string inputs from the user and echos back the first 
*  letter of each wored until a word starting with 9 is encountered.
*/

#include <stdio.h>

#define MAX_LENGTH 25

main()
{
	char word[MAX_LENGTH];

	printf("\nEnter one word at a time followed by <enter>\n");	
	printf("> ");
	scanf("%s", word);

	while (word[0] != '9')
	{
		printf("\n%s starts with the letter %c\n", word, word[0]);
		printf("> ");
		scanf("%s", word);	
	}
}



