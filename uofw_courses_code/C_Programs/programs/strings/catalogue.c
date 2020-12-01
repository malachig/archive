/* Author: Malachi Griffith
*  Date: Nov. 17 2002  
*  Purpose: Seperates a product ID of the form, ATL1203S14, into its components
*  City: ATL   Product: 1203   Qualifiers: S14
*/

#include <stdio.h>
#include <string.h>
#define MAX_LENGTH 20
main()
{
	char id_value[MAX_LENGTH];
	char city[10];
	char product[10];
	char size[5];

	int i = 0;
	int start, end;

	/* Get the product ID number */
	printf("Enter the product ID > ");
	scanf("%s", id_value);	

	while (id_value[i] >= 'A' && id_value[i] <= 'Z')
		i++;

	start = i;

	strncpy(city, id_value, i);
	city[i] = '\0';	
	printf("City: ");
	puts(city);	
	
	while (id_value[i] >= '0' && id_value[i] <= '9')
		i++;
	
	/* Use & ADDRESS OPERATOR TO COPY FROM A PARTICULAR STARTING POSITION */
	strncpy(product, &id_value[start], i - start);
	product[i] = '\0';
	printf("Product: ");
	puts(product);
	
	start = i;

	/* Rest of string is size value. */
	strcpy(size, &id_value[start]);
	printf("Size: ");
	puts(size);
 

}



