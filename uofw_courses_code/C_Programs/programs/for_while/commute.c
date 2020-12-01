/* commute.c */

/* Author: Malachi Griffith
*  Date: Oct. 25 2002
*  Purpose: Calcuates amount of subsidy for peop who carpool.
*	Min passenger efficiency (in passenger-kilometers per liter) = 25
*	P = (n * s) / l
*  	n = number of passengers.
*	s = distance travelled in km.
*	l = litres of gas used.
*  Subsidy is $0.08 per passenger-km where passenger-km = (n * s)
*/

#include <stdio.h>
#define REBATE 0.08
main()
{
	int n;
	double s;
	double l;
	double pass_km;
	double P;
	double subsidy = 0.00;

	FILE *input_data;
	input_data = fopen("commute.dat", "r");

	printf("\nPassengers\tWeekly Commute\tGasoline\tEfficiency\tWeekly");
	printf("\n\t\t(km)\t\tConsumption(L)\t(pass km / L)\tSubsidy($)");

	fscanf(input_data, "%d%lf%lf", &n, &s, &l);
	
	while(n != 0)
	{
		pass_km = n * s;
		P = pass_km / l;
		
		if (P > 25)
			subsidy = pass_km * REBATE;
		printf("\n%d\t\t\t%.0f\t%.1f\t%.1f\t%.2f",
			 n, s, l, P, subsidy);
		
		fscanf(input_data, "%d%lf%lf", &n, &s, &l);
	}
	printf("\n\n");	
	fclose(input_data);
}



