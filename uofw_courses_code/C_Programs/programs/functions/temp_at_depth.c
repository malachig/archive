/* Author: Malachi Griffith
*  Date: Oct. 24 2002
*  Purpose: Calculates the temperature below the surface of the earth
*  (distance in km) in both farenheit and celsius.
*/

#include <stdio.h>

double temp_c(double distance);
double temp_f(double temp_c);

main()
{
	double depth;
	double Tc, Tf;

	printf("\nEnter the distance below the surface in km > ");
	scanf("%lf", &depth);

	Tc = temp_c(depth);
	Tf = temp_f(Tc);

	printf("\n\nThe temperature is %.2f degrees C  or %.2f degrees F\n",
		Tc, Tf);
}

/*
*  Function: temp_c
*/
double
temp_c(double distance)
{
	double T;
	T = (10 * distance) + 20;
	return(T);
}

/* 
*  Function: temp_f
*/
double
temp_f(double temp_c)
{
	double T;
	
	T = (1.8 * temp_c) + 32;
	return(T);
}
 

