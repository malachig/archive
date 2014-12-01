/* Author: Malachi Griffith
*  Date: Dec. 7 2002  
*  Purpose: Allocation of arrays with calloc()
*/

#include <stdio.h>
#include <stdlib.h>

typedef struct{
	char name[10];
	double diameter;
	int moons;
	double orbit_time;
	double rotation_time;
	} planet_t;

/* Function Prototype */
int scan_planet(planet_t *plnp);

int
main()
{
	char *string1;
	int *array_of_nums;	
	planet_t *array_of_planets;

	int str_siz, num_nums, num_planets, i;

	printf("Enter string length and string > ");
	scanf("%d", &str_siz);

	string1 = (char *)calloc(str_siz, sizeof(char));
	scanf("%s", string1);

	printf("\nHow many numbers? > ");
	scanf("%d", &num_nums);

	array_of_nums = (int *)calloc(num_nums, sizeof(int));
	array_of_nums[0] = 5;

	for(i = 1; i < num_nums; ++i)
		array_of_nums[i] = array_of_nums[i-1]*i;

	printf("\nEnter number of planets and planet data > ");
	scanf("%d", &num_planets);
	array_of_planets = (planet_t *)calloc(num_planets, sizeof(planet_t));

	for (i = 0; i < num_planets; i++)
		scan_planet(&array_of_planets[i]);
}

int 
scan_planet(planet_t *plnp)  /* output - address of planet_t structure
			      * to fill. */
{
	int result;

	result = scanf("%s%lf%d%lf%lf", plnp->name,
					&plnp->diameter,
					&plnp->moons,
					&plnp->orbit_time,
					&plnp->rotation_time);

	if (result == 5)
		result = 1;
	else if (result != EOF)
		result = 0;

	return (result);
}

