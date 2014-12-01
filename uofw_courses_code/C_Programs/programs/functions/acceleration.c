/* Author: Malachi Griffith
*  Date: Oct. 24 2002 
*  Purpose: Calculate acceleration using the formula,
*		a = (v_final - v_initial) / t
*/

#include <stdio.h>

/* Function prototype */
double acceleration(double initial_velocity, double final_velocity,
		    double time);
double time_to_stop(double initial_velocity, double acceleration);

main()
{
	double Vi,
	       Vf,
	        t,
	        a;
	double stopping_time;

	printf("\nEnter the initial velocity > ");
	scanf("%lf", &Vi);
	printf("Enter the final velocity > ");
	scanf("%lf", &Vf);
	printf("Enter the time for this change to occur > ");
	scanf("%lf", &t);

	a = acceleration(Vi, Vf, t);
	stopping_time = time_to_stop(Vi, a);

	printf("\nThe acceleration is %.2f m/s^2", a);
	printf("\nThe time to stop is %.2f seconds\n\n", stopping_time);
	
}


/*
*  Function: acceleration
*/
double
acceleration(double initial_velocity, double final_velocity, double time)
{
	double accel;

	accel = (final_velocity - initial_velocity) / time;

	return(accel);
}

/*
*  Function: time_to_stop
*/
double
time_to_stop(double initial_velocity, double acceleration)
{
	double final_v = 0;
	double time;

	time = (final_v - initial_velocity) / acceleration;

	return(time);
}
	
	

