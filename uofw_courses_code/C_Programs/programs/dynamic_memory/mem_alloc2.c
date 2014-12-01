/* Author: Malachi Griffith
*  Date: Nov. 21 2002
*  Purpose: To illustrate the value of using memory allocation this example
*  will declare only a pointer for each neccessary possible record.  As 
*  each record is needed, then the space will be allocated by dynamic 
*  memory allocation.  
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>  /* Needed for dynamic memory allocation */

#define NAMESIZE 51
#define NUM_MPS 10 

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

/* Declare an array of pointers to the type Employee, redefining the typedef */
typedef Employee * Employee_ptr;

FILE *infile;

/* Function Prototypes */
void input(Employee_ptr m_o_p[], int *size);
void display(Employee_ptr m_o_p[], int size);

main()
{
	int size;
	Employee_ptr m_o_p[NUM_MPS];
	input(m_o_p, &size);
	display(m_o_p, size);
}

/*
*  Function: input()
*/
void
input(Employee_ptr m_o_p[], int *size)
{
	int i = 0;
	Employee temp;
	
	infile = fopen("mop.dat", "r");
	
	while(!feof(infile))
	{
		/* Dynamic Memory Allocation */
	  	m_o_p[i] = (Employee_ptr)malloc(sizeof(Employee));

		fscanf(infile, "%s%s%s%d%lf", temp.emp_name.lastname,
			temp.emp_name.firstname,
			temp.emp_name.middlename,
			&temp.employee_id,
			&temp.salary);
		*m_o_p[i] = temp;
		i++;
	}
	*size = i - 1; 
	fclose(infile);
}

void display (Employee_ptr m_o_p[], int size)
{
	int i;
	for (i = 0; i < size; i++)
	{	
		printf("Name = %s %s %s\n", ( *m_o_p[i]).emp_name.firstname,
			(*m_o_p[i]).emp_name.middlename, 
			(*m_o_p[i]).emp_name.lastname);

		printf("ID = %d\n", (*m_o_p[i]).employee_id);
		printf("Salary = $%.2lf\n", (*m_o_p[i]).salary);
		printf("\n\n");
	}
}
