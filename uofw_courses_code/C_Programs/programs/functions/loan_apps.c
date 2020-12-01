/* Author: Malachi Griffith
*  Date: Oct. 14 2002
*  Purpose: A program for approving or rejecting loan applications
*/

#include <stdio.h>

#define APPROVE 1
#define REJECT 0
#define PREFERRED 1
#define LOW_BARRIER 10
#define HIGH_BARRIER 12

/* Function prototype declarations */

void process_app(int[], int, int, double, double);
void display_message(int);
int decide_loan(int, int, double, double);

main()
{
	int term;			/* Term of loan */
	int status;			/* Customer status */
	int counters[] = {0, 0};	/* Approval/rejection counters */

	double income;			/* Customer yearly income */
	double amount;			/* Amount of loan */

/* Read data from file and process application */

	while(!feof(stdin))
	{
	scanf("%lf%d%lf%d\n", &income, &term, &amount, &status);
	process_app(counters, term, status, income,amount);
	}

/* Display total number of approvals and rejections */

	printf("Total number of approvals: %4d\n", counters[APPROVE]);
	printf("Total number of rejections: %4d\n", counters[REJECT]);
}

/* This function processes the application by calling 
*  decide_loan, displays the appropriate message and 
*  increments the counters.
*/

void process_app(int counters[], int term, int status, double income,
	         double amount)
{
	int approval;
	
	approval = decide_loan(term, status, income, amount);
	
	display_message(approval);

	if (approval)
		counters[APPROVE]++;
	else
		counters[REJECT]++;
	return;
}

/* This function decides whether to approve a loan application
*  or not.
*/

int decide_loan(int term, int status, double income, double amount)
{
	int decision;

	if ((amount / term) <= (income / HIGH_BARRIER))
		decision = APPROVE;
	else
		if ((status == PREFERRED) &&
	 	   ((amount / term) <= (income / LOW_BARRIER)))
		   decision = APPROVE;
		else
		   decision = REJECT;

	return (decision);
}

/* This function displays the appropriate message */

void display_message(int decision)
{
	if (decision)
		printf("The application is approved.\n");
	else
		printf("The application is rejected.\n");

	return;
}
	

