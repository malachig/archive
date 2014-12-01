/* Assignment 4: Question 1 (15 Marks) sorting_structures.c */

/* Author: Malachi Griffith
*  Date: Nov. 24 2002 
*  Purpose: This program will use a data structure to store information
*  from a data file.  The data will consist of 1.) A Property Class 
*  2.) A Selling Price.  and 3.) The Commission Earned.
*  Once entered, the data will be sorted by property class (ascending) and 
*  within property class, by selling price (ascending).
*  Finally the sorted list of data will be printed out to the file,
*  "sorting_structures.print".
*/

#include <stdio.h>

#define MAX_PROPERTIES 25  /* Max # of entries in data file */

#define COM_RATE1 0.045    /* Commission rate for class 1 property */
#define COM_RATE2 0.050    /* Commission rate for class 2 property */
#define COM_RATE3 0.060    /* Commission rate for class 3 property */

/* Define a structure to store property sale record */
typedef struct {
		int property_class;        /* A value of 1, 2, or 3. */
		double selling_price;      /* Dollar value of sale */
		double commission_earned;  /* Based on comission rate for
					    * that class times sale price */
	}Sale;

/* Function Prototypes */
void read_array(Sale sales_record[], int *size);
void sort_array(Sale sales_record[], int size);
int find_min (Sale sales_record[], int start, int end);
void display_array(Sale sales_record[], int size);

main()
{
	int size_local; /* Number of properties actually found in datafile */

	Sale sales_record[MAX_PROPERTIES]; /* Array of structures */
	
	/* Call user defined functions */
	read_array(sales_record, &size_local);
	sort_array(sales_record, size_local);
	display_array(sales_record, size_local);
}

/*
*  Function: read_array()
*  Reads data from the data file "assign4.dat" and enters them into an
*  array of structures of type Sale.  
*  Pre: The sales_record array was defined in main.  
*  Post: The array is filled as long as their is data in the datafile
*  	 and the number of entries found is counted and updated in main
*  	 as size_local. 
*/
void
read_array(Sale sales_record[], int *size)
{
	FILE *input_data;  /* File pointer for input datafile */
	int i = 0; 	   /* Counter in while loop and loop control 
			    * variable in for loop.
	
	/* Open the file for reading data */
	input_data = fopen("assign4.dat", "r");

	/* Get an initial priming set of data from the file. */
	fscanf(input_data, "%d%lf", &sales_record[i].property_class,
	       &sales_record[i].selling_price);

	/* Continue to get data and count the entries until the file ends. */
	while(!feof(input_data))
	{
		i++;
		fscanf(input_data, "%d%lf", &sales_record[i].property_class,
	               &sales_record[i].selling_price);
	}
	*size = i;  /* Tells main() how many entries need to be processed */

	/* Determine and input the commission earned for each record */
	/* Commission Rate applied depends on the property class */
	for (i = 0; i < *size; i++)
	{
		if (sales_record[i].property_class == 1)
			sales_record[i].commission_earned = 
			(sales_record[i].selling_price * COM_RATE1);
	
		else if (sales_record[i].property_class == 2)
			sales_record[i].commission_earned =
			(sales_record[i].selling_price * COM_RATE2);
		
		else if (sales_record[i].property_class == 3)
			sales_record[i].commission_earned =
			(sales_record[i].selling_price * COM_RATE3);

		else 
			printf("\n**Error in Data File - Invalid Class**\n");
	}
	fclose(input_data);
}

/*
*  Function: sort_array()
*  Sorts the array in ascending order accoring to property class AND selling
*  price.
*  Pre: The array of structures, 'sales_record' and number of entries, 'size'
*  	are defined. 
*  Post: Rearranges the order of the structures in the array until they 
*	 are completely sorted.
*/
void 
sort_array(Sale sales_record[], int size)
{
	Sale temp_record;  /* Temp structure to facilitate swapping */
	int index; 	   /* Loop control variable */
	int current; 	   /* Position of current min found by find_min() */

	/* Examine every position in the array and at each position place 
	*  the current min value at that position unless there already */
	for (index = 0; index < size -1; index++)
	{
		current = find_min(sales_record, index, size);
		
		/* As long as the current min differs from the position
		*  being considered, do a swap. */
		if (index != current)
		{
			temp_record = sales_record[index];
			sales_record[index] = sales_record[current];
			sales_record[current] = temp_record;
		}
	}
}

/*
*  Function: find_min()
*  This function is a helper function that finds the current minimum value
*  from a given starting point sent to the end of the array.
*  Pre: The sales_record array, starting point, and size of the array are
*       defined.
*  Post: The position of the current minimum value within the subarray is
* 	 returned to sort_array as an integer.
*/
int 
find_min (Sale sales_record[], int start, int end)
{
	int current;  /* Current minimum position */
	int i; 	      /* Loop control variable */
 
	current = start; /* Initially set current min to the first position
			 *  in the subarray being considered. */

	/* Scan through the subarray */
	for (i = start + 1; i < end; i++)
	{
	/* If the property class is LOWER than the current minimum then 
	*  this position is definitely a new minimum position */ 
	   if (sales_record[i].property_class < 
	       sales_record[current].property_class)
			current = i;

	/* If the property class is the SAME as the current minimum then
	*  this position is only a new minimum if the selling price is 
	*  also lower!. */
	   else if ((sales_record[i].property_class ==
	   	    sales_record[current].property_class) && 
		    (sales_record[i].selling_price <
		     sales_record[current].selling_price))
			current = i;
	}

	return (current);  /* Return the new current minimum to sort_array */
}

/*
*  Function: display_array()
*  Simply prints out a summary to the output file "sorting_structures.print"
*  Pre:  The sales_record array and number of entries, size are defined.
*  Post:  Prints out summary to output file. 
*/
void
display_array(Sale sales_record[], int size)
{
	int i; 			/* Loop control variable */
	FILE *output_data;	/* File pointer variable */ 

	/* Open the output file for writing */
	output_data = fopen("sorting_structures.print", "w");
	
	/* Print the table headings. */
	fprintf(output_data, "Property Class\t\tSelling Price ($)");
	fprintf(output_data, "\tCommission Earned ($)\n");

	/* Print out the data itself */
	for(i = 0; i < size; i++)
		fprintf(output_data, "%7d\t\t\t%13.2f\t\t%13.2f\n",
			sales_record[i].property_class,
			sales_record[i].selling_price,
			sales_record[i].commission_earned);

	fprintf(output_data, "\n\n");

	fclose(output_data);
}
