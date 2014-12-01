/* Author: Malachi Griffith
*  Date: Nov. 23
*  Purpose: Gets one line of data from standard input.  Returns an empty
*  string on end of file.  If data line will not fit in allocated space,
*  stores portion that does fit and discards rest of input line.
*/

#include <stdio.h>
#include <string.h>

char * scanline(char *dest, int dest_len);

main()
{
	char string1[81];
	char string2[81];
	char ch;
	
	scanline(string1, 40);
 
	puts(string1);





}

char *
scanline(char *dest,	/* output - destination string */
	 int dest_len)	/* input - space available in dest */
{
	int i, ch;

	/* Gets next line one character at a time. */
	i = 0;

	for (ch = getchar();
	     ch != '\n' && ch != EOF && i < dest_len - 1;
	     ch = getchar())

		dest[i++] = ch;

	dest[i] = '\0';

	/* Discards any characters that remain on input line */
	while(ch != '\n' && ch != EOF)
		ch = getchar();

	return(dest);
}


