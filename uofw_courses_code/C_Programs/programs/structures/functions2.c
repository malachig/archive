/* Author: Malachi Griffith
*  Date: Nov. 20 2002
*  Purpose: Using functions along with structures.  Note that an array
*  of structures can be easily read in.  HOWEVER in this case it is not
*  an array.  You therefore need to use a pointer to a structure.
*/

#include <stdio.h>
#include <string.h>

#define NAMESIZE 51

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
void input(Employee *m_o_p);
void display(Employee m_o_p);

main()
{
	Employee m_o_p;
	
	input(&m_o_p);

	display(m_o_p);
}

/*
*  Function: input()
*/
void
input(Employee *m_o_p)
{
	infile = fopen("mop.dat", "r");	
	
	fscanf(infile, "%s%s%s%d%lf", (*m_o_p).emp_name.lastname,
		(*m_o_p).emp_name.firstname,
		(*m_o_p).emp_name.middlename,
		&(*m_o_p).employee_id,
		&(*m_o_p).salary);
	fclose(infile);
}

/*
* Function: display()
*/
void 
display(Employee m_o_p)
{
		printf("Name = %s %s %s\n",
			m_o_p.emp_name.firstname,
			m_o_p.emp_name.middlename,
			m_o_p.emp_name.lastname);

		printf("ID = %d\n", m_o_p.employee_id);
		printf("Salary is $%.2f\n", m_o_p.salary);
		printf("\n\n");
}




