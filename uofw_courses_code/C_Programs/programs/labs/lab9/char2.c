/* Author: Malachi Griffith
*  Date: Oct. 31 2002
*  Purpose:  Dealing with character values.
*/

#include <stdio.h>

main()
{
	char ch;
	int i = 0;
	int j;

	char letters[50] = {' '};

	FILE *input_data;
	input_data = fopen("char.dat", "r");

	while(!feof(input_data))
	{
		fscanf(input_data, "%c", &ch); 
		if (ch == ' ')
			continue;
	
		if (ch >= 97 && ch <= 122)
			ch -= 32;
		letters[i] = ch;
		i++;
	
		while (ch != '\n')
		{
			fscanf(input_data, "%c", &ch);
			if (ch >= 65 && ch <= 90)
				ch += 32;	
			letters[i] = ch;
			i++;	
		}
	}
	fclose(input_data);

	/* PRINT OUT THE CHARACTERS IN THE ARRAY */
	printf("\n");
	for(j = 0; j <= i; j++)
		printf("%c", letters[j]);

}



