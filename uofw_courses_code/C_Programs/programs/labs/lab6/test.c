/* Author: Malachi Griffith
*  Date: Oct. 10 2002
*  Purpose: Calculate a persons net pay.
*	    Gross pay is based on hours worked and wage plus any overtime
*	    they work.  This will be regular pay times 1.5 for hours over 40.
*	    They must also pay 15% tax on their total gross pay.
*/

#include <stdio.h>

/* Function Prototypes */
double gross_pay(double pay_rate, double hours);
double tax(double gross_pay);
double net_pay(double gross_pay, double tax);

main()
{
	double pay_rate;
	double hours;
	double gross;
	double taxes;
	double net;

	printf("\n\nEnter the number of hours you worked > ");
	scanf("%lf", &hours);
	printf("Enter the pay rate > ");
	scanf("%lf", &pay_rate);

	gross = gross_pay(pay_rate, hours);

	taxes = tax(gross);

	net = net_pay(gross, taxes);

	return(0);
}

/* Function gross pay */

double gross_pay(double pay_rate, double hours)
{
	double pay;
	double overtime;
	double overtime_pay;
	double total_pay;

	if (hours < 0)
		printf("Not a valid entry for hours");
	else if (hours <= 40)
		pay = hours * pay_rate;
	else
		{
		pay = 40 * pay_rate;
		overtime = (hours - 40);
		overtime_pay = overtime * (pay_rate * 1.5);
		total_pay = pay + overtime_pay;
                }

	return(total_pay);
}
		
/* Function tax */

double tax(double gross)
{
	double total_tax;

	total_tax = gross * 0.15;

	return(total_tax);
}

/* Function net pay */

double net_pay(double gross_pay, double tax)
{
	double grand_total;

	grand_total = gross_pay - tax;

	return(grand_total);
}
 
