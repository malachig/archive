/* Author: Malachi Griffith
*  Date: Oct. 22 2002
*  Purpose: Converts race time into speed in ft/s and m/s.
*/

#include <stdio.h>
#define FEET_PER_MILE 5280
#define FEET_PER_KILOMETER 3282


/* Function prototypes */
void input_data(int *minutes, double *seconds);
double feet_per_second(int minutes, double seconds);
double meters_per_second(double fps);

main()
{
	int m;
	double s;	
	double feet, meters;

	input_data(&m, &s);
	feet = feet_per_second(m, s);
	meters = meters_per_second(feet);

	printf("\nYour speed was: %f fps or %f mps\n\n", feet, meters); 
}

/*
*  Function: input_data
*/
void
input_data(int *minutes, double *seconds)
{
	int min;
	double sec;
	
	printf("\nEnter the number of minutes you took > ");
	scanf("%d", &min);
	printf("Enter the number of seconds you took > ");
	scanf("%lf", &sec);

	*minutes = min;
	*seconds = sec;
}


/* 
* Function feet_per_second
*/

double
feet_per_second(int minutes, double seconds)
{
	int distance;
	double feet_per_second;
	double time_in_s;

	distance = FEET_PER_MILE;
	time_in_s = seconds + (minutes * 60);
	feet_per_second = (double) distance / time_in_s;

	return(feet_per_second);
}

/*
*  Function: meters_per_second
*/
double 
meters_per_second(double fps)
{
	double feet_per_meter;
	double meters_per_second;

	feet_per_meter = FEET_PER_KILOMETER / 1000;
	meters_per_second = fps / feet_per_meter;

	return(meters_per_second);
}


	















