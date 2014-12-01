/* Author: Malachi Griffith
*  Date: Oct. 21 2002
*  Purpose: Computes duration of a rockets flight and its height above
*  ground when it reaches its target.
*/

#include <stdio.h>
#include <math.h>

#define G 32.17  /* Gravitational constant */
#define DEGREES_IN_A_RADIAN 57.295

/* Function Prototypes */
void input_data (double *theta, double *distance, double *velocity);
void calculations (double theta, double distance, double velociy,
		   double *time, double *height);
void print_out (double time, double height);

main()
{
	double theta_value;
	double dist;
	double v;
	double t;
	double h;
	
	input_data(&theta_value, &dist, &v);
	calculations(theta_value, dist, v, &t, &h);
	print_out(t, h);




}

/* Function: input_data */
void
input_data(double *theta, double *distance, double *velocity)
{
	double theta_local;
	double distance_local;
	double velocity_local;

	printf("\nEnter the theta value in degrees > ");
	scanf("%lf", &theta_local);
	*theta = theta_local;

	printf("\nEnter the distance value > ");
	scanf("%lf", &distance_local);
	*distance = distance_local;

	printf("\nEnter the velocity > ");	
	scanf("%lf", &velocity_local);
	*velocity = velocity_local;
	printf("\n\n");

}
	 
/* Function Calculations */
void
calculations (double theta, double distance, double velocity,
	      double *time, double *height)
{
	double time_local;
	double height_local;

/*	theta *= DEGREES_IN_A_RADIAN;*/
	time_local = (distance) / (velocity * (cos(theta)));
	*time = time_local;

	height_local = (velocity * (sin(theta)) * time_local) - 
			((G * time_local * time_local) / 2);
	*height = height_local;
}

/* Function print_out */
void 
print_out (double time, double height)
{
	printf("The time value is: %.2f\n", time);
	printf("The height value is: %.2f\n\n", height);
}
