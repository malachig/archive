/* Author: Malachi Griffith
*  Date: Nov. 10 2002
*  Purpose: Scans sales figures for one year and stores them in a table 
*  organized by salesperson and quarter.  Displays the table and the annual
*  totals for each person and the totals for each quarter.
*/

#include <stdio.h>

#define SALES_FILE "sales1.dat"  /* name of sales data file */
#define NUM_QUARTERS 4
#define  NUM_SALES_PEOPLE 5

typedef enum {fall, winter, spring, summer} quarter_t;

/* Function Prototypes */
int scan_table(double sales[] [NUM_QUARTERS], int num_rows);
void sum_rows(double col_sum[], double sales[] [NUM_QUARTERS], int num_rows);
void sum_columns(double col_sum[], double sales[] [NUM_QUARTERS],
		 int num_rows);
void display_table(double sales[] [NUM_QUARTERS], 
		   const double person_totals[], 
		   const double quarter_totals[],
		   int num_rows);
void initialize(double sales[] [NUM_QUARTERS], int num_rows, double value);
void display_quarter(quarter_t quarter);

main()
{
	double sales [NUM_SALES_PEOPLE] [NUM_QUARTERS];  /* Table of sales */
	double person_totals[NUM_SALES_PEOPLE];	 	 /* Row totals */
	double quarter_totals [NUM_QUARTERS];		 /* Column totals */
	int status;

	status = scan_table(sales, NUM_SALES_PEOPLE);
	if (status == 1)
	{
	
		sum_rows(person_totals, sales, NUM_SALES_PEOPLE); 
		sum_columns(quarter_totals, sales, NUM_SALES_PEOPLE);
	
		display_table(sales, person_totals, quarter_totals,
			      NUM_SALES_PEOPLE);

	}
	return(0);
}


/*
*  Scans the sales data from SALES_FILE and computes and stores the sales
*  results in the sales table.  Flags out-of-range data and data format errors
*  Post:  Each entry of sales represents the sales total for a particular
*	  salesperson and quarter.
*	  Returns 1 for a successful table scan, 0 for error in scan.
*  Calls: Initialize to initialize table to all zeros.
*/
int
scan_table (double sales[] [NUM_QUARTERS], /* output */
	    int num_rows)		   /* input */
{
	double trans_amt;	/* transaction amount */
	int trans_person;	/* salesperson number */
	quarter_t quarter;	/* sales quarter */
	FILE *sales_filep;	/* file pointer to sales file */
	int valid_table = 1;	/* data valid so far */
	int status;		/* input status */
	char ch;		/* one character in bad line */

	/* Initialize table to all zeros */
	initialize(sales, num_rows, 0.0);
	
	/* Scan and store the valid sales data */
	sales_filep = fopen(SALES_FILE, "r");

	for(status = fscanf(sales_filep, "%d%d%lf", &trans_person,
			    &quarter, &trans_amt);
	    status == 3 && valid_table;
	    status = fscanf(sales_filep, "%d%d%lf", &trans_person,
			    &quarter, &trans_amt))
	{
		if (fall <= quarter && quarter <= summer &&
		    trans_person >= 0 && trans_person < num_rows)
		{
			sales[trans_person][quarter] += trans_amt;
		}
		else
		{
			printf("Invalid salesperson or quarter -- \n");
			printf("  person is %d, quarter is ", trans_person);
			
			display_quarter(quarter);
			printf("\n\n");
		
			valid_table = 0;
		}
	}

	if (!valid_table)  /* error already processed */
	{
		status = 0;
	}
	else if (status == EOF)  /* end of data file without error */
	{
		status = 1;
	}
	else
	{
		printf("Error ins sales data format.  Revise data.\n");
		printf("ERROR HERE >>> ");
		for (status = fscanf(sales_filep, "%c", &ch);
		     status == 1 && ch != '\n';
		     status = fscanf(sales_filep, "%c", &ch))
			
			printf("%c", ch);
		
		printf(" <<<\n");
		status = 0;
	}
	return (status);
}


/*
*  Stores value in all elements of sales.
*  Pre:  Value is defined and num_rows is the number of rows in sales.
*  Post: All elements of sales have desired value.
*/
void
initialize (double sales[] [NUM_QUARTERS],	/* output */
	    int num_rows,			/* input */
	    double value)			/* input */
{
	int row;
	quarter_t quarter;

	for(row = 0; row < num_rows; ++row)
		for(quarter = fall; quarter <= summer; ++quarter)
			sales[row] [quarter] = value;
}

/*
*  Displays the sales table data in table form along with row and column
*  sums.
*  Pre: sales , person_totals, quarter_totals, and num_rows are defined.
*  Post:  Values stored in the three arrays are displayed.
*/

void
display_table(double sales[] [NUM_QUARTERS],	/* input */
	      const double person_totals[],	/* input */
	      const double quarter_totals[],	/* input */
	      int num_rows)
{
	int person;
	quarter_t quarter;

	/* Display heading */
	printf("%34cSALES Summary\n%34c----- -------\n\n", ' ', ' ');
	printf("%12s%5c", "Salesperson", ' ');

	for (quarter = fall; quarter <= summer; ++quarter)
	{
		display_quarter(quarter);
		printf("%8c", ' ');
	}
	printf("TOTAL\n");
	printf("--------------------------------------");
	printf("--------------------------------------\n");

	/* Display Table */
	for (person = 0; person < num_rows; ++person)
	{
		printf("%6d%4c", person, ' ');
		for (quarter = fall; quarter <= summer; ++quarter)
			printf("%6c%8.2f", ' ', sales[person][quarter]);
		printf("%6c%8.2f\n", ' ', person_totals[person]);	
	}
	
	printf("--------------------------------------");
	printf("--------------------------------------\n");
	printf("QUARTER TOTALS ");
	for (quarter = fall; quarter <= summer; ++quarter)
		printf("%9.2f%5c", quarter_totals[quarter], ' ');
	printf("\n");
}

/*
*  Display an enumeration constant of type quarter_t
*/

void
display_quarter(quarter_t quarter)
{
	switch(quarter)
	{
		case fall:	printf("Fall");
				break;

		case winter:	printf("Winter");
				break;
		
		case spring:	printf("Spring");
				break;

		case summer:	printf("Summer");
				break;

		default:	printf("Invalid quarter %d", quarter);
	}
}

/*
*  Function: sum_rows.  Calculates the total sales amount for each person
*/

void 
sum_rows(double col_sum[], double sales[] [NUM_QUARTERS], int num_rows)
{

	int i, j;  
	
	for(i = 0; i < num_rows; i++)
		for(j = 0; j < NUM_QUARTERS; j++)
		{
			col_sum[i] += sales[i][j];
		}

}

/*
*  Function: sum_columns.  Calculates the total sales amount for each
*  quarter of the year.
*/

void sum_columns(double col_sum[], double sales[] [NUM_QUARTERS],
		 int num_rows)
{

	int i, j;

	for (j = 0; j < NUM_QUARTERS; j++)
		for (i = 0; i < num_rows; i++)
		{
			col_sum[j] += sales[i][j];
		}

}

