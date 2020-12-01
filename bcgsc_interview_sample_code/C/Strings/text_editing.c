/* Author: Malachi Griffith
*  Date: Nov. 23 2002 
*  Purpose: Performs text editing operations on a source string 
*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define MAX_LEN 100
#define NOT_FOUND -1

/* Function Prototypes */
char *delete(char *source, int index, int n);
char *do_edit(char *source, char command);
char get_command(void);
char *insert(char *source, const char *to_insert, int index);
int pos(const char *source, const char *to_find);

int
main(void)
{
	char source[MAX_LEN];
	char command;

	printf("Enter the source string:\n");
	gets(source);

	for(command = get_command();
	    command != 'Q';
	    command = get_command())
	{
		do_edit(source, command);
		printf("New source: %s\n\n", source);
	}
	
	printf("String after editing: %s\n", source);
	return(0);
}

/*
*  Returns source after deleting n characters beginning with source[index].
*  If source is too short for full deletion, as many characters are deleted
*  as possible.
*  Pre: All parameters are defined and strlen(source) - index - n < MAX_LEN
*  Post:  source is modified and returned
*/
char *
delete(char *source,	/* input/output - string from which to delete part */
       int index,	/* input - index of first char to delete */
       int n)		/* input - number of chars to delete */
{
	char rest_str[MAX_LEN];	/* Copy of souce substring following 
				*  characters to delete */

	/* If there are no characters in source following portion to 
	delete, delete rest of string */

	if (strlen(source) <= index + n)
	{
		source[index] = '\0';
	}

	/* Otherwise, copy the portion following the portion to delete
	*  and place it in source beginning at the index position */
	else
	{
		strcpy(rest_str, &source[index + n]);
		strcpy(&source[index], rest_str);
	}
	
	return(source);
}

/*
*  Performs the edit operation specified by command
*  Pre: command and source are defined.
*  Post:  After scanning additional information needed, performs a 
*   	  deletion (command = 'D') or insertion (command = 'I') or
*	  finds a substring ('F') and displays result; returns
*	  (possibly modified) source.
*/
char *
do_edit(char *source,  /* input/output - string to modify or search */
	char command)  /* input - character indicating operation */
{
	char str[MAX_LEN];	/* work string */
	int index;

	switch (command)
	{
	case 'D':
		printf("String to delete> ");
		gets(str);
		index = pos(source, str);
		if (index == NOT_FOUND)
			printf("'%s' not found\n", str);
		else
			delete(source, index, strlen(str));
		break;

	case 'I':
		printf("String to insert > ");
		gets(str);
		printf("Position of insertion > ");
		scanf("%d", &index);
		insert(source, str, index);
		break;

	case 'F':
		printf("String to find > ");
		gets(str);
		index = pos(source, str);
		if (index == NOT_FOUND)
			printf("'%s' not found\n", str);
		else
			printf("'%s' found at position %d\n", str, index);
		break;

	default:
		printf("Invalid edit command '%c'\n", command);
	}

	return(source);
}

/*
*  Prompt for and get a character representing an edit command and 
*  convert it to uppercase.  Return the uppercase character and ignore 
*  rest of input line.
*/
char
get_command(void)
{
	char command, ignore;

	printf("Enter D (delete), I(insert), F(Find), or Q (Quit) > ");
	scanf(" %c", &command);

	do
		ignore = getchar();
	while (ignore != '\n');

	return (toupper(command));
}

/*
*  Returns source after inserting to_insert at postion index of source.
*  If source[index] doesn't exist, adds to_insert at end of source.
*  Pre: all parameters are defined, space available for source is 
*  enough to accomodate insertion, and strlen(source) - index - n < MAX_LEN
*  Post: source is modified and returned
*/
char *
insert(char *source, 		/* input/output - target of insertion */
       const char *to_insert,	/* input - string to insert */
       int index)		/* input - position where to_insert is 
				*  is to be inserted*/
{
	char rest_str[MAX_LEN]; /* copy of rest of source beginning 
				* with source[index] */

	if (strlen(source) <= index)
	{
		strcat(source, to_insert);
	}
	else
	{
		strcpy(rest_str, &source[index]);
		strcpy(&source[index], to_insert);
		strcat(source, rest_str);
	}

	return(source);
}

/*
*  Returns index of first occurence of to_find in source or 
*  value of NOT_FOUND if to_find is not in source.
*  Pre: both parameters are defined.
*/
int 
pos(const char *source,  /* input -  string in which to look for to_find */
    const char *to_find)  /* input - string to find. */
{
	int i = 0, find_len, found = 0, position;
	char substring[MAX_LEN];

	find_len = strlen(to_find);

	while (!found && i <= strlen(source) - find_len)
	{
		strncpy(substring, &source[i], find_len);
		substring[find_len] = '\0';

		if (strcmp(substring, to_find) == 0)
			found = 1;
		else
			++i;
	}

	if (found)
		position = i;
	else 
		position = NOT_FOUND;

	return(position);
}
