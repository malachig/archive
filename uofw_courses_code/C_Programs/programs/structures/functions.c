/* Author: Malachi Griffith
*  Date: Nov. 20 2002
*  Purpose: Using functions along with structures.  Note that an array
*  of structures can be easily read in.
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

	for (i = 0; i < NUM_MPS; i++)
	{
		fscanf(infile, "%s%s%s%d%lf", m_o_p[i].emp_name.lastname,
			m_o_p[i].emp_name.firstname,
			m_o_p[i].emp_name.middlename,
			&m_o_p[i].employee_id,
			&m_o_p[i].salary);
	}
	fclose(infile);
}

/*
* Function: display()
*/
void 
display(Employee m_o_p[])
{
	int i;

	for (i = 0; i < NUM_MPS; i++)
	{
		printf("Name = %s %s %s\n",
			m_o_p[i].emp_name.firstname,
			m_o_p[i].emp_name.middlename,
			m_o_p[i].emp_name.lastname);

		printf("ID = %d\n", m_o_p[i].employee_id);
		printf("Salary is $%.2f\n", m_o_p[i].salary);
		printf("\n\n");
	}
}




