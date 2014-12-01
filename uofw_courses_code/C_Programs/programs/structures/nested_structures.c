/* Author: Malachi Griffith
*  Date: Nov. 20 2002
*  Purpose: Demonstrate the use of structures nested within other 
*  structures.
*/

#include <stdio.h>
#include <string.h>
#define NAMESIZE 51

typedef struct {
		char first[NAMESIZE];
		char middle[NAMESIZE];
		char last[NAMESIZE];
		} Name;

typedef struct {
		Name emp_name;  /* name */
		int employee_id;       /* ID */
		double salary;	       /* salary */
	       } Employee;

main()
{
	Employee prime_minister;

	Employee leader_opposition;

	char line[NAMESIZE];

	prime_minister.salary = 216000.00;
	strcpy(prime_minister.emp_name.first, "Jean");
	strcpy(prime_minister.emp_name.middle, "JJ");
	strcpy(prime_minister.emp_name.last, "Cretien");


	printf("%s makes $%.2f a year.\n", prime_minister.emp_name.last,
		prime_minister.salary);

	leader_opposition.salary = 131000.00;
	leader_opposition.employee_id = 2;
	gets(line);

	strcpy(leader_opposition.emp_name.last, line);
	
	printf("Will %s ever get those $%.3f???\n", 
		leader_opposition.emp_name.last,
	        prime_minister.salary);
}
