/* compare_strings.c */

/* Author: Malachi Griffith
*  Date: Nov 21 2002
*  Purpose: Use a function called compare() that will accept tow strings as
*  its two arguments. Assume that their max length is 30.  If the two strings
*  are different it will return 0, if they are the same it will return 1.
*  NO LIBRARY FUNCTIONS ARE ALLOWED FOR THIS PROGRAM 
*/

#include <stdio.h>
#define MAX_SIZE 30

/* Function Prototypes */
int compare(char str1[], char str2[], int size1, int size2);

main()
{
	char string1[MAX_SIZE] = {0};
	char string2[MAX_SIZE] = {0};
	
	char ch;
	int i = 0;

	int result;
	int size1, size2;

	printf("\nPlease enter the first string >\n");

	scanf("%c", &ch);
	
	while (ch != '\n')
	{
		string1[i] = ch;
		i++;
		scanf("%c", &ch);
	} 
	size1 = i;

	i = 0;
	printf("\nPlease enter the second string >\n");

	scanf("%c", &ch);

	while (ch != '\n')
	{
		string2[i] = ch;
		i++;
		scanf("%c", &ch);
	}
	size2 = i;

	result = compare(string1, string2, size1, size2);

	if (result == 0)
		printf("The strings are not the same\n");
	else if (result == 1)
		printf("The strings are the same\n");

}

int
compare(char str1[], char str2[], int size1, int size2)
{
	int i;
	int compare_size;
	int test = 1;

	if (size1 != size2)
		{
		test = 0;
		return(test);
		}

	compare_size = size1;

	for (i = 0; i <= compare_size; i++)
	{
		if (str1[i] != str2[i])
		{
			test = 0;
			break;	
		}
	}

	return(test);
}
		







