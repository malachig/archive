/* Author: Malachi Griffith
*  Date: Nov. 20 2002
*  Purpose: Illustrate the basic usage of the 'to' operator in combination
*  with nested structures.
*/

#include <stdio.h>
#define NAMESIZE 51

typedef struct{
	char firstname[NAMESIZE];
	char lastname[NAMESIZE];
	char middlename[NAMESIZE];
	}  Name;

typedef struct{
		Name emp_name;
		int employee_id;
		double salary;
		} Employee;

main()
{
	Employee prime_minister = {"Jean", "Cretien", "JJ", 1, 262000.00};

	Employee *pm_ptr;

	pm_ptr = &prime_minister;

	printf("Salary = $%.2f\n", pm_ptr->salary);
	printf("Name = %s %s %s\n", pm_ptr->emp_name.firstname,
		pm_ptr->emp_name.middlename, pm_ptr->emp_name.lastname);
}
