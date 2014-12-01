/* Author: Malachi Griffith
*  Date: Nov. 21 2002
*  Purpose: To illustrate the value of using memory allocation this example
*  will declare a static sized table.  The next example mem_alloc2.c will 
*  use memory allocation to do the same thing.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>  /* Needed for dynamic memory allocation */

#define NAMESIZE 51
#define NUM_MPS 3

typedef struct {
		char firstname[NAMESIZE];
		char lastname[NAMESIZE];
		char middlename[NAMESIZE];
	} Name;

typedef struct{
		Name emp_name;
		int employee_id;
		double salary;
	} Employee;

FILE *infile;

/* Function Prototypes */
void input(Employee m_o_p[]);
void display(Employee m_o_p[]);

main()
{
	Employee m_o_p[NUM_MPS];
	input(m_o_p);
	display(m_o_p);
}

/*
*  Function: input()
*/
void
input(Employee m_o_p[])
{
	int i;
	
	infile = fopen("mop.dat", "r");
	
	for(i = 0; i < NUM_MPS; i++)
	{
		fscanf(infile, "%s%s%s%d%lf", m_o_p[i].emp_name.lastname,
			m_o_p[i].emp_name.firstname,
			m_o_p[i].emp_name.middlename,
			&m_o_p[i].employee_id,
			&m_o_p[i].salary);
	}
	fclose(infile);
}

void display (Employee m_o_p[])
{
	int i;
	for (i = 0; i < NUM_MPS; i++)
	{	
		printf("Name = %s %s %s\n", m_o_p[i].emp_name.firstname,
			m_o_p[i].emp_name.middlename, 
			m_o_p[i].emp_name.lastname);

		printf("ID = %d\n", m_o_p[i].employee_id);
		printf("Salary = $%.2lf\n", m_o_p[i].salary);
		printf("\n\n");
	}
}
