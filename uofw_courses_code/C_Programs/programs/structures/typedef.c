/* Author: Malachi Griffith
*  Date: Nov. 20 2002
*  Purpose: Demonstrate the basic use of structures
*/

#include <stdio.h>
#include <string.h>
#define NAMESIZE 51

typedef struct {
		char name[NAMESIZE];  /* name */
		int employee_id;       /* ID */
		double salary;	       /* salary */
	       } Employee;

main()
{
	Employee prime_minister = {"Jean Cretien", 1, 262000.0};

	Employee leader_opposition;

	char line[NAMESIZE];

	printf("%s makes $%.2f a year.\n", prime_minister.name,
		prime_minister.salary);

	leader_opposition.salary = 131000.00;
	leader_opposition.employee_id = 2;
	gets(line);

	strcpy(leader_opposition.name, line);
	
	printf("Will %s ever get those $%.3f???\n", leader_opposition.name,
	       prime_minister.salary);
}
