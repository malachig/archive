/* Author: Malachi Griffith
*  Date: Nov. 20 2002
*  Purpose: Define an array of structures to hold the records for
*  each category of employee.  Assume that the following would read 
*  a datafile which contains, each on a single line, the names, employee
*  id number, and salaries of the Members of Parliament.
*/

#include <stdio.h>
#include <string.h>
#define NAMESIZE 51
#define NUM_MPS 2

typedef struct{
		char firstname[NAMESIZE];
		char lastname[NAMESIZE];
		char middlename[NAMESIZE];
		} Name;

typedef struct {
		Name emp_name;
		int employee_id;
		double salary;
		} Employee;

main()
{
	Employee m_o_p[NUM_MPS];

	int i;

	for(i = 0; i < NUM_MPS; i++)
	{
		gets(m_o_p[i].emp_name.firstname);
		gets(m_o_p[i].emp_name.lastname);
		gets(m_o_p[i].emp_name.middlename);
		scanf("%d", &m_o_p[i].employee_id);
		scanf("%lf\n", &m_o_p[i].salary);
	}

	for (i = 0; i < NUM_MPS; i++)
	{
		printf("Name = %s %s %s \n",
				m_o_p[i].emp_name.firstname,
				m_o_p[i].emp_name.middlename,
				m_o_p[i].emp_name.lastname);
		
		printf("ID = %d\n", m_o_p[i].employee_id);
		printf("Salary = $%.2lf\n", m_o_p[i].salary);
		printf("\n\n");
	}
}
