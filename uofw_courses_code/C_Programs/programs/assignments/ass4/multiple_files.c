/* Assignment 4: Question 2 (15 Marks) multiple_files.c */

/* Author: Malachi Griffith
*  Date: Nov. 24 2002 
*  Purpose: This program will use a data structure to store information
*  from a data file.  The data will consist of 1.) A Property Class 
*  2.) A Selling Price.  and 3.) The Commission Earned.
*  Once entered, the data will be printed to output files depending on 
*  their class (i.e., the sales for class 1 will be printed to "out1.dat",
*  class 2 entries to "out2.dat" and class 3 to "out3.dat". 
*  Input data will be read from assign4.dat.  Finally, a function will 
*  be used to read the output files, display the contents and summarize
*  the total commission earned for the class.  The user will enter the 
*  desired class file to display. 
*/

#include <stdio.h>
#include <string.h>

#define MAX_PROPERTIES 25   /* Max number of entries in data file */
#define MAX_FILENAME 20     /* Max length of filename specified by user */

#define COM_RATE1 0.045   /* Commission rate for class 1 property */
#define COM_RATE2 0.050   /* Commission rate for class 2 property */
#define COM_RATE3 0.060   /* Commission rate for class 3 property */

/* Define a struture to store property_sale records */
typedef struct {
		int property_class;	  /* Class of property (1,2,3)*/
		double selling_price;     /* Sale price in dollars */
		double commission_earned; /* Commision Rate * Sale Price */
	}Sale;


/* Function Prototypes */
void make_files(Sale sales_record[]);
void display_files(char filename[]);

main()
{
	Sale sales_record[MAX_PROPERTIES];  /* Array or structures */
	char filename[MAX_FILENAME]; 	    /* Array of char to store file
					    *  name entered by user */
	char answer[5] = {"Y"};		    /* Array of char to store answer
					    *  from user about repeating the
					    *  program or not */

	/* Call user defined function to create three output files */	
	make_files(sales_record);

	/* Istructions for user displayed to screen */	
	printf("\nOutput files have been named out1.dat for class 1, ");
	printf("out2.dat for class 2, etc.");

	/* Continue showing files until the user is finished */
	while (!strcmp(answer, "Y") || !strcmp(answer, "y"))
	{
	  printf("\n\nEnter the name of the file to display > ");
	  gets(filename);
	  
	  /* If the filename does not match one of the three, repeat loop
	  *  and ask them to enter it again!. */
	  if ((strcmp(filename, "out1.dat") && strcmp(filename, "out2.dat")
	      && strcmp(filename, "out3.dat")))
		{
		printf("** Invalid Filename! **\n");
		continue;
		}
	
	  /* Call user defined function to display file summary to screen */
	  display_files(filename);
	
	  /* Ask the user if they wish to continue or exit */	
	  printf("Would you like to view another file? (Y or N) > ");
	  gets(answer); 
	}
}

/*
*  Function: make_files()
*  This function simply reads the data from "assign4.dat" into an array of
*  structures, calculates the commission for each record, and then goes 
*  through each record and prints it to one of three output files,
*  "out1.dat", "out2.dat", or "out3.dat" depending on its class (1, 2, or 3).
*  Pre: The array of structures 'sales_record' is defined.
*  Post: Three output datafiles are generated.
*/
void
make_files(Sale sales_record[])
{
	/* Define File pointer variables for the input file and each 
	*  of the output files. */ 
	FILE *input_data;
	FILE *out1;
	FILE *out2;
	FILE *out3;

	int i = 0;  /* Counter for number of entries */
	int size;   /* Number of entries actually found */
	
	/* Open the data files for reading or writing */
	input_data= fopen("assign4.dat", "r");
	out1 = fopen("out1.dat", "w");
	out2 = fopen("out2.dat", "w");
	out3 = fopen("out3.dat", "w");

	/* Get an initial set of data from the input file */
	fscanf(input_data, "%d%lf", &sales_record[i].property_class,
	       &sales_record[i].selling_price);

	/* Continue to get data until the file ends */
	while(!feof(input_data))
	{
	  i++;
	  fscanf(input_data, "%d%lf", &sales_record[i].property_class,
	         &sales_record[i].selling_price);
	}
	size = i;  /* Number of data records */

	/* Determine and input the commission earned for each record */
	for (i = 0; i < size; i++)
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

		else    /* Display error message to standard output */
			printf("\n**Error in Data File - Invalid Class**\n");
	}

	/* Direct each record to its corresponding output file */
	/* Class 1 records go to out1.dat, etc. Each file is specified
	 * by its file pointer variable */
	for (i = 0; i < size; i++)
	{
	  if (sales_record[i].property_class == 1)
		fprintf(out1, "%.2f\t%.2f\n", 
			sales_record[i].selling_price,
			sales_record[i].commission_earned);
	  else if (sales_record[i].property_class == 2)
		fprintf(out2, "%.2f\t%.2f\n", 
			sales_record[i].selling_price,
			sales_record[i].commission_earned);
	  else if (sales_record[i].property_class == 3)
		fprintf(out3, "%.2f\t%.2f\n", 
			sales_record[i].selling_price,
			sales_record[i].commission_earned);
	}
	/* Close all files after done using them */
	fclose(input_data);
	fclose(out1);
	fclose(out2);
	fclose(out3);
}

/*
*  Function: display_files()
*  Simply displays a summary of an output file to the screen in tabular 
*  format with simple headings.
*  Pre: Filename specifed by user is sent as an argument from main().
*  Post: Screen display.
*/
void
display_files(char filename[])
{
	double sale_price, 	/* Local variable for selling price */
	       commission; 	/* Local variable for commission earned */ 

	double total_commission = 0; /* Accumulator for commission earned
				     *  by the current property class */
	FILE *input_data;	/* File pointer variable */

	/* Open the file specified by the user for reading */ 
	input_data = fopen(filename, "r");	

	/* Print out title and headings to screen */	
	printf("\nSummary for Class Stored in '%s'\n\n", filename);
	printf("Selling Price ($)\tCommission Earned ($)\n");

	/* Get a priming set of data */
	fscanf(input_data, "%lf%lf", &sale_price, &commission);

	/* Display the data and get another set until EOF */
	while(!feof(input_data))
	{
		printf("%13.2f\t\t%13.2f\n", sale_price, commission);
		total_commission += commission;	
		fscanf(input_data, "%lf%lf", &sale_price, &commission);
	}
	printf("\nThe total commission earned for this class was $%.2f\n\n",
		total_commission);

	fclose(input_data);
}
