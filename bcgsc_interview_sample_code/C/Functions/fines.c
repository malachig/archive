/* Assignment 2: Question 3 Parts A and B (15 Marks)  fines.c */

/* Author: Malachi Griffith
*  Date: Oct. 13 2002
*  Purpose: Calculates fines for parking violations.  The program will
*  read data from an input file called "fines.dat" which consists of 
*  the number of each type of violation incurred by one person.  Using
*  the conditions provided in the question, the total value of the fines
*  will be added and printed out.  Four functions will be called by the
*  main function to accomplish this task:
*	1.) base_fine   - finds the base fine value for each person.
*	2.) extra_fines - tallies the additional fines incurred.
*	3.) total_fine  - calculates the grand total.
*	4.) print_out - prints out the results.	 
*/

#include <stdio.h>

#define COST_A 10.00  	/* Fine for each Type A violation. */
#define COST_B 20.00	/* Fine for each Type B violation. */
#define COST_C 30.00 	/* Fine for each Type C violation. */
#define COST_D 45.00 	/* Fine for each Type D violation. */


/* Function Prototypes */
double base_fine  (int type_a, int type_b, int type_c, int type_d);
double extra_fine (int type_a, int type_b, int type_c, int type_d);
double total_fine (double base_fine, double extra_fine);
void print_out    (double base_fine, double extra_fine, double total_fine,
	           int person_count);

main()
{
	int type_a,	/* Number of Type A violations per person */
	    type_b,	/* Number of Type B violations per person */
	    type_c, 	/* Number of Type C violations per person */
	    type_d; 	/* Number of Type D violations per person */

	int person_count = 0;  /* The current person being processed */

	double base; 	/* Total base fine for that person */
	double extra;	/* Total additional fines for that person */	
	double total;   /* Grand total fines for that person */

	/* Define and open the data file */
	FILE * input_fines;
	input_fines = fopen("fines.dat", "r");
	
	fscanf(input_fines, "%d%d%d%d", &type_a, &type_b, &type_c, &type_d);
	while (!feof(input_fines))
	{
		person_count++;	

		/* Calls to User Defined Functions */	
		base = base_fine(type_a, type_b, type_c, type_d);

		extra = extra_fine(type_a, type_b, type_c, type_d);

		total = total_fine(base, extra);

		print_out(base, extra, total, person_count);
		
		/* retrieve data for the next person from the data file */
		fscanf(input_fines, "%d%d%d%d", &type_a, &type_b, &type_c,
		       &type_d);
	}

	printf("\n\n");
	fclose(input_fines);
}


/* Function base_fine */
/* Calculates the base fine in $, based on the number of each violation.
*  Pre:  type_a, type_b, type_c, and type_d are defined
*  Post: A value for the total base fine is returned to the main function,
*        as base_fine.
*/

double
base_fine(int type_a, int type_b, int type_c, int type_d)
{
	/* Local Variables */

	double 	fine_a, 	/* Total value of Type A fines */ 
	   	fine_b,		/* Total value of Type B fines */
		fine_c,		/* Total value of Type C fines */
		fine_d;		/* Total value of Type D fines */
	double base_fine; 	/* Total of all base fines added together */

	/* Multiply the cost for each violation by the number of them */
	fine_a = type_a * COST_A;
	fine_b = type_b * COST_B;
	fine_c = type_c * COST_C;
	fine_d = type_d * COST_D;
	
	base_fine = fine_a + fine_b + fine_c + fine_d;
	
	return(base_fine);
}


/* Function extra_fine */
/* Calculates the extra fines according to the conditional statements below.
*  Pre: type_a, type_b, type_c, type_d are defined.
*  Post: A value for the total extra_fine is returned to the main function.
*/

double
extra_fine(int type_a, int type_b, int type_c, int type_d)
{
	/* Local Variables */

	int extra_d;		/* Number of type D fines exceeding 3 */
	int total_abc;		/* Number of Type A, B, and C added up */
	double extra_fine = 0;  /* Total cost of all extra_fines incurred */

	/* Extra penalty if number of Type A, B or C fines exceeds 10 */
	if (type_a > 10)
		extra_fine += 50.00;
	if (type_b > 10)
		extra_fine += 50.00;
	if (type_c > 10)
		extra_fine += 50.00;

	/* Extra penalty for Each Type D fine over 3 */
	if (type_d > 3)
		{
		extra_d = (type_d - 3) * 20.00;
		extra_fine += extra_d;
		}
	
	/* Extra penalty if none of Type A, B, or C exceeds 10 individually
	*  but the sum of these three types exceeds 20 */
	total_abc = type_a + type_b + type_c;
	
	if (total_abc > 20 && type_a <= 10 && type_b <= 10 && type_c <=10)
		extra_fine += 75.00;

	return(extra_fine);
}



/* Function total_fine */
/* Calculates the total fines incurred by each person 
*  Pre: base_fine and extra_fine are defined
*  Post: A value for the total_fine incurred is returned to the main
*  function.
*/

double
total_fine(double base_fine, double extra_fine)
{
	double total_fine; 	/* Grand total value of fines per person */

	total_fine = base_fine + extra_fine;
	
	return(total_fine);
}


/* Function print_out */
/* Prints the results to the screen or output file in a readable manner.
*  Pre: base_fine, extra_fine, and total_fine are defined.
*  No values are returned to the main function
*/

void
print_out(double base_fine, double extra_fine, double total_fine,
          int person_count)
{
	printf("\nThe base fine incurred by person %d is:    $%.2f", 
	       person_count, base_fine);
	printf("\nThe extra fines incurred by person %d are: $%.2f",
	       person_count, extra_fine);
	printf("\nThe total fines incurred by person %d are: $%.2f\n",
	       person_count, total_fine); 
}

