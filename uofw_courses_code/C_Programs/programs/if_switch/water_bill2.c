/* water_bill2.c */

/* Author: Malachi Griffith
*  Date: Sept. 28 2002
*  Purpose: Computes and prints a water bill given an unpaid balance and previous and
*  current meter readings.  Bill includes a demand charge of $35.00, a use
*  charge of $1.10 per thousand gallons, and a surcharge of $2.00 if there is 
*  an unpaid balance.
*/

#include <stdio.h>

#define DEMAND_CHG 35.00
#define PER_1000_CHG 1.10
#define LATE_CHG 2.00
#define OVERUSE_CHG_RATE 2.0  /* Double use charge as non-conservation penalty */
#define CONSERV_RATE 95  /* Percent of last year's use allowed this year */

/* Function Prototypes */
void instruct_water(void);

double comp_use_charge (int previous, int current, int use_last_year);

double comp_late_charge(double unpaid);

void display_bill(double late_charge, double bill, double unpaid);

int
main()
{
	int 	previous; /* input - previous meter reading in thoudands of gallons */
	int 	current;  /* input - meter reading from current quarter */
	int	use_last_year;  /* Use for same quarter last year */
	double 	unpaid;   /* input - unpaid balance */
	double  bill;	  /* output - water bill */
	int	used;	  /* Thousands of gallons used this quarter */
	double	use_charge;  /* Charge for actual water use */
	double	late_charge; /* Charge for non-payment of part of previous bill */

	/* Display user instructions */
	
	instruct_water();
	
	/* Get data: unpaid balance, previous and current meter readings. */

	printf("Enter unpaid balance > $");
	scanf("%lf", &unpaid);
	printf("Enter previous meter reading > ");
	scanf("%d", &previous);
	printf("Enter the current meter reading > ");
	scanf("%d", &current);
	printf("Enter your usage for this quarter last year > ");
	scanf("%d", &use_last_year);

	/* Compute user charge */
	use_charge = comp_use_charge(previous, current, use_last_year);

	/* Determine applicable late charge */
	late_charge = comp_late_charge(unpaid);
	
	/* Figure bill */
	bill = DEMAND_CHG + use_charge + unpaid + late_charge;

	/* Print bill. */
	display_bill(late_charge, bill, unpaid);
	
	return(0);
}


/*
*  Displays user instructions 
*/

void
instruct_water(void)
{
	printf("This program figures a water bill ");
	printf("based on the demand charge\n");
	printf("($%.2f) and a $%.2f per 1000 ", DEMAND_CHG, PER_1000_CHG);
	printf("gallons use charge.\n\n");
	printf("A $%.2f surcharge is added to ", LATE_CHG);
	printf("accounts with an unpaid balance.\n");
	printf("\nEnter unpaid balance, previous ");
	printf("and current meter readings\n");
	printf("on seperate lines after the prompts.\n");
	printf("Press <return> or <enter> after ");
	printf("typing each number.\n");
	printf("You must also enter your water use in the ");
	printf("quarter last year\n\n");
}

/*
* Computes use charge
* Pre: previous, current, and use_last_year are defined.
*/

double
comp_use_charge(int previous, int current, int use_last_year)
{
	int used;  /* Gallons of water used in thousands */
	double use_charge;  /* charge for actual water use */

	used = current - previous;

	if (used <= CONSERV_RATE / 100.0 * use_last_year)
	{
		/* Conservation guidlines met */
		use_charge = used * PER_1000_CHG;
	}
	else
	{
		printf("Use charge is at %.2f times ", OVERUSE_CHG_RATE);
		printf("normal rate since use of\n");
		printf("%d units exceeds %d percent ", used, CONSERV_RATE);
		printf("of last year's %d-unit use.\n", use_last_year);
		use_charge = used * OVERUSE_CHG_RATE * PER_1000_CHG;
	}

	return (use_charge);
}


/*
*  Computes late charge.
*  Pre: unpaid is defined.
*/

double
comp_late_charge(double unpaid)
{
	double late_charge;  /* charge for non-payment of part of previous balance */

	if (unpaid > 0)
		late_charge = LATE_CHG; /* Assess late charge on unpaid balance */
	else
		late_charge = 0.0;

	return (late_charge);
}


/*
*  Displays late charge if any and bill.
*  Pre: late_charge, bill, and unpaid are defined.
*/
void
display_bill(double late_charge, double bill, double unpaid)
{
	if (late_charge > 0.0)
	{
		printf("\nBill includes $%.2f late charge", late_charge);
		printf(" on unpaid balance of $%.2f\n", unpaid);
	}
	printf("\nTotal due = $%.2f\n", bill);
} 

