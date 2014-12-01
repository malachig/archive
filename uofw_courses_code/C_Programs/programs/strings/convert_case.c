/* Author: Malachi Griffith
*  Date: Nov. 23 2002  
*  Purpose: Converts the lowercase letters of its string argument to 
*  uppercase leaving other characters unchanged.
*/

#include <stdio.h>
#include <ctype.h>

#define STRSIZ 80

/* Function Prototype */
char * string_toupper(char *str);
int string_greater(const char *str1, const char *str2);

main()
{






}

/*
* Function: string_toupper()
*/
char *
string_toupper(char *str)
{
	int i;

	for (i = 0; i < strlen(str); ++i)
		if (islower(str[i]))
			str[i] = toupper(str[i]);

	return(str);
}

/* 
*  Function: string_greater
*  Compares two strings of up to STRSIZ characters ignoring the case of
*  the letters.  Returns the value 1 if str1 should follow str2 in an
*  alphabetized list; otherwise returns 0
*/
int 
string_greater (const char *str1, const char *str2)
{
	char s1[STRSIZ], s2[STRSIZ];

	/* Copies str1 and str2 so string_toupper can modify copies */
	strcpy(s1, str1);
	strcpy(s2, str2);

	return(strcmp(string_toupper(s1), string_toupper(s2)) > 0);
}

