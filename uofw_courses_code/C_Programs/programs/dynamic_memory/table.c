/* Author: Malachi Griffith
*  Date: Nov. 21 2002 
*  Purpose: Illustrate the use of tables constructed as an array of 
*  structures.  
*/

#include <stdio.h>
#include <stdlib.h>  /* Needed for dynamic memory allocation */
#include <string.h>

#define MAX_CHAR 8
#define TABLE_SIZE 5

typedef struct {
		char part_number[MAX_CHAR];
		int quantity_on_hand;
		double unit_price;
	} Inventory;

/* Function Prototypes */
void Read_table(Inventory []);
void Print_table(Inventory []);
void Sort_table(Inventory []);

main()
{
	Inventory inventory_table[TABLE_SIZE];
	char key[MAX_CHAR];
	
	Read_table(inventory_table);
	printf("Before sorting:\n");
	Print_table(inventory_table);
	Sort_table(inventory_table);
	printf("After sorting:\n");
	Print_table(inventory_table);
/*
	while (!(feof(stdin)))
	{
		gets(key);
		printf("Searching for the key %s\n", key);
		Search_table(inventory_table, key);
		scanf("\n");
	}
*/
}

void
Read_table(Inventory inventory_table[TABLE_SIZE])
{
	int row;
	int ichar;

	/* Read in table values */
	
	for(row = 0; row < TABLE_SIZE; row++)
	{
		for (ichar = 0; row < MAX_CHAR -1; ichar++)
			scanf("%c", &inventory_table[row].part_number[ichar]);
	
		inventory_table[row].part_number[MAX_CHAR - 1] = '\0';

		scanf("%d", &inventory_table[row].quantity_on_hand);
		scanf("%lf", &inventory_table[row].unit_price);
		scanf("\n");
	}
}

void
Print_table(Inventory inventory_table[TABLE_SIZE])
{
	int row;

	for (row = 0; row < TABLE_SIZE; row++)
	{
		printf("%s\t%d\t%.2f\n",
			inventory_table[row].part_number,
			inventory_table[row].quantity_on_hand,
			inventory_table[row].unit_price);
	}
	printf("\n");
}

/*
*  Function: Sort_table()
*  This function performs a bubble sort on a table of type Inventory 
*/
void
Sort_table(Inventory inventory_table[TABLE_SIZE])
{
	Inventory temp;

	int limit;
	int exchange;
	int index;
	int flag;

	index = 0;
	limit = TABLE_SIZE - 1;
	flag = 1;
	exchange = 0;

	while (flag)  /* ie. while flag is 1 not 0. */
	{
		flag = 0; /* Will be changed back to 1 if exchange occurs */

		for (index = 0; index < limit; index++)
		{
			/* Compare first part number to next one, etc. */
			if(strcmp(inventory_table[index].part_number,
			   inventory_table[index + 1].part_number) > 0)
			   {
			   /* If bigger it 'floats' up (swaps with next) */ 
			   temp = inventory_table[index];
			   inventory_table[index]=inventory_table[index+1];
			   inventory_table[index+1] = temp;
			   flag = 1; /* Because an exchange has occured */
			   /* Update position of last exchange */ 
			   exchange = index; 
			   }
		}
		/* Position of last exchange for that pass.  If an exchange 
		*  occured, another pass will take place, etc. until no
		*  exchanges occur. */
		limit = exchange;
	}
}
