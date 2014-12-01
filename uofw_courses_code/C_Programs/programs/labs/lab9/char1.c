/* Author: Malachi Griffith
*  Date: Oct. 31 2002
*  Purpose:  Dealing with character values.
*/

#include <stdio.h>

main()
{
	char ch;

	FILE *input_data;
	input_data = fopen("char.dat", "r");

	while(!feof(input_data))
	{
		fscanf(input_data, "%c", &ch); 
		if (ch == ' ')
			continue;
	
		if (ch >= 97 && ch <= 122)
			ch -= 32;
		printf("%c", ch);
	
		while (ch != '\n')
		{
			fscanf(input_data, "%c", &ch);
			if (ch >= 65 && ch <= 90)
				ch += 32;	
			printf("%c", ch);
		}

	printf("\n");
	}
	fclose(input_data);
}



