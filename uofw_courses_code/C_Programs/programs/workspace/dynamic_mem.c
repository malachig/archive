#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NAMESIZE 51
#define NUM_MPS 10

typedef struct{
	char firstname[NAMESIZE];
	char lastname[NAMESIZE];
	char middle[NAMESIZE];
		}Name;

typedef struct{
	Name emp_name;
	int employee_id;
	double salary;
	} Employee;

typedef Employee *EmployeePtr;

FILE *infile;

/* Function Prototypes */
void input(EmployeePtr m_o_p[], int *size);
void display(EmployeePtr m_o_p[], int size);

main()
{
	int size;
	EmployeePtr m_o_p[NUM_MPS];

	input(m_o_p, &size);
	display(m_o_p, size);
}

void
input (EmployeePtr m_o_p[], int *size)
{
	int i = 0;
	Employee temp;
	infile = fopen("mop.dat", "r");

	while(!feof(infile))
	{
	  m_o_p[i] = (EmployeePtr)malloc(sizeof(Employee));

	  fscanf(infile, "%s%s%s%d%lf\n",
	  	 temp.emp_name.lastname,
		 temp.emp_name.firstname,
		 temp.emp_name.middle,
		 &temp.employee_id,
		 &temp.salary);

	  *m_o_p[i] = temp;
	  i++;
	}
	*size = i;
	fclose(infile);
}

void
display (EmployeePtr m_o_p[], int size)
{
	int i;

	for (i = 0; i < size; i++)
	{
		printf("Name = %s %s %s\n",
		       m_o_p[i]->emp_name.firstname,
		       m_o_p[i]->emp_name.middle,
		       m_o_p[i]->emp_name.lastname);
		printf("ID = %d\n",
		       m_o_p[i]->employee_id);
		printf("Salary = $%.2f",
		       m_o_p[i]->salary);
	printf("\n\n");
	}
}
	  
