/* char_replace.c */

/* Author: Malachi Griffith
*  Date: Nov. 7 2002 
*  Purpose: Take a string from the user at standard input, then ask the user
*  for a character, search for this character and replace it with a 
*  a new character specified by the user.
*/

#include <stdio.h>

#define MAX_SIZE 81

/* Function prototypes */
int find_char(char string[], char ch);
void replace_character(char string[], char replacement, int position);

main()
{
	char string[MAX_SIZE];
	char search_char;
	char replace_char;
	int position;

	printf("Please enter a sentence > \n");
	gets(string);

	printf("Please enter a character to search for > ");
	scanf("%c", &search_char);	

	printf("Please enter a replacement for this character > ");
	scanf(" %c", &replace_char);

	position = find_char(string, search_char); 	

	printf("\nThe character is at position %d", position + 1);

	replace_character(string, replace_char, position);

	printf("\nThe edited sentence is as follows:\n");
	puts(string);	
	printf("\n");
}


/* Function find_char */
int 
find_char(char string[], char ch)
{
	int i;
	int position;

	for (i = 0; i < MAX_SIZE; i++)
	{
		if (string[i] == ch)
			{
			position = i;
			break;
			}
	}

	return(position);
}	

/* Function replace_character */
void 
replace_character(char string[], char replacement, int position)
{
	string[position] = replacement;	
}





