/* Author: Malachi Griffith
*  Date: Dec 3 2002 
*  Purpose: Converts measurements given in one unit to any other unit of
*  the same category that is listed in the database file, units.dat.
*  Handles both names and abbreviations of units.
*/

#include <stdio.h>
#include <string.h>

#define NAME_LEN 30  	/* Storage allocated for a unit name */
#define ABBREV_LEN 15	/* Storage allocated for a unit abbreviation */
#define CLASS_LEN 20	/* Storage allocated for a measurement class */
#define NOT_FOUND -1	/* Value indicating unit not found */
#define MAX_UNITS 20	/* Maximum number of different units handled */

typedef struct {		/* Unit of measurement type */	
	char name[NAME_LEN]; 	/* Character string such as "milligrams" */
	char abbrev[ABBREV_LEN]; /* Shorter character string such as "mg" */
	char class[CLASS_LEN];	/* Character string such as "pressure",
				*  "distance", "mass", etc.  */
	double standard;	/* Number of standard units equivalent to
				*  this unit. */
	} unit_t;

int fscan_unit(FILE *filep, unit_t *unitp);
void load_units(int unit_max, unit_t units[], int *unit_sizep);
int search(const unit_t units[], const char *target, int n);
double convert(double quantity, double old_stand, double new_stand);

int
main(void)
{
	unit_t units[MAX_UNITS];  /* Units classes and conversion factors */
	int num_units;		  /* number of elements of units in use */
	char old_units[NAME_LEN], /* units to convert (name or abbrev) */
	     new_units[NAME_LEN]; /* units to convert to (name or ...) */
	int status;		  /* input status */
	double quantity;	  /* value to convert */
	int old_index;		  /* index of units element where old_units
				   * found. */
	int new_index;		  /* index where new_units found */

	/* Load units of measurement database */
	load_units(MAX_UNITS, units, &num_units);

	/* Convert quantities to desired units until data format error 
	*  including error code returned when q is entered to quit */
	printf("Enter a conversion problem or q to quit.\n");
	printf("To convert 25 kilometers to miles, you would enter\n");
	printf("> 25 kilomteres miles\n");
	printf("  or, alternatively, \n");
	printf("> 25 km mi\n> ");

	for (status = scanf("%lf%s%s", &quantity, old_units, new_units);
	     status == 3;
	     status = scanf("%lf%s%s", &quantity, old_units, new_units))
	{
		printf("Attempting aconversion of %.4f %s to %s . . .\n",
		       quantity, old_units, new_units);

		old_index = search(units, old_units, num_units);
		new_index = search(units, new_units, num_units);

		if (old_index == NOT_FOUND)
			printf("Unit %s not in database\n", old_units);
		else if (new_index == NOT_FOUND)
			printf("Unit %s not in database\n", new_units);
		else if (strcmp(units[old_index].class,
				units[new_index].class) != 0)
			printf("Cannot convert %s (%s) to %s (%s)\n",
				old_units, units[old_index].class,
				new_units, units[new_index].class);
		else 
			printf("%.4f%s = %.4f %s\n", quantity, old_units,
			       convert(quantity, units[old_index].standard,
			       units[new_index].standard),
			       new_units);

		printf("\nEnter a conversion problem or q to quit.\n> ");
	}
	return(0);
}

/*
*  Gets data from a file to fill output argument
*  Returns standard error code:  1 => successful input, 0 => error,
* 				 negative EOF value => end of file
*/
int 
fscan_unit(FILE *filep,	  /* input - file pointer */
	   unit_t *unitp) /* output - unit_t structure to fill */
{
	int status;

	status = fscanf(filep, "%s%s%s%lf", unitp->name, 
					    unitp->abbrev,
					    unitp->class,
					    &unitp->standard);

	if (status == 4)
		status = 1;
	else if (status != EOF)
		status = 0;

	return(status);
}

/*
*  Opens database file units.dat and gets data to place in units until
*  end of file is encountered.  Stops input prematurely if there are more
*  than unit_max data values in the file or if invalid data is encountered.
*/
void
load_units (int	unit_max,	/* input - declared size of units */
	    unit_t units[],	/* output - array of data */
	    int *unit_sizep)	/* output - number of data values stored
				*  in units. */
{
	FILE *inp;
	unit_t data;
	int i, status;

	/* Gets database of units from file */
	inp = fopen("units.dat", "r");
	i = 0;

	for (status = fscan_unit(inp, &data);
	     status == 1 && i < unit_max;
	     status = fscan_unit(inp, &data))
	{
		units[i++] = data;
	}
	fclose(inp);

	/* Issue error message on premature exit */
	if (status == 0)
	{
		printf("\n*** Error in data format ***\n");
		printf("*** Using first %d data values ***\n", i);
	}
	else if (status != EOF)
	{
		printf("\n***Error: too much data in file ***\n");
		printf("*** Using first %d data values ***\n", i);
	}

	/* Send back size of used portion of array */
	*unit_sizep = i;
}

/*
*  Searches for target key in name and abbrev components of first n
*  elements of array units
*  Returns index of structures containing target or NOT_FOUND
*/
int
search(const unit_t units[],  /* array of unit_t structures to search */
       const char *target,    /* key searched for in name and abbrev
			       * components */
       int n)		      /* number of array elements to search */
{
	int i,
	found = 0,	/* whether or not target has been found */
	where;	   	/* index where target found or NOT_FOUND */

	/* Compare name and abbrev components of each element to target */
	i = 0;

	while (!found && i < n)
	{
		if (strcmp(units[i].name, target) == 0 ||
		    strcmp(units[i].abbrev, target) == 0)
			found = 1;
		else
			++i;
	}
	
	/* Return index of element containing target or NOT_FOUND */
	if (found)
		where = i;
	else 
		where = NOT_FOUND;

	return(where);
}

/*
*  Converts one measurement to another given the representation of both 
*  in a standard unit.  For example, to convert 24 feet to yards given a
*  standard unit of inches: quantity = 24, old_stand = 12 (there are 12
*  inches in a foot), new_stand = 36 (there are 36 inches in a yard),
*  results is 24 * 12 / 36 which equals 8
*/
double
convert (double quantity,	/* value to convert */
	 double old_stand,	/* number of standard units in one of 
				*  quantity's original units */
	 double new_stand)	/* number of standard units in 1 new unit */
{
	return (quantity * old_stand / new_stand);
}
	







