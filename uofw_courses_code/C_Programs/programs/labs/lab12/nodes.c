/* nodes.c */

/* Author: Malachi Griffith
*  Date: Nov. 27 2002 
*  Purpose: A batch program to calculate the retirement benefits for its
*  employees for each month.  Only full time employees can have retirement
*  benefits.  The benefit is 5% of the pay earned during the month.  Casual
*  employees are not entitled to retirement benefits.
*/

#include <stdio.h>
#include <string.h>

typedef struct node{
			char emp_id[10];
			int emp_status;
			double pay_rate;
			double hours_worked;
			double benefit;
			struct node *link;
		} Node;

typedef Node *NodePtr;  

/* Function Prototypes */
void PrintList(NodePtr node_ptr); 

main()
{
	FILE *input_data;

	NodePtr temp;
	NodePtr head = NULL;
	NodePtr last = NULL;

	char temp_id[10];
	int temp_status;
	double temp_pay;
	double temp_hours;
	double temp_benefit;

	int i;

	input_data = fopen("nodes.dat", "r");

	for (i = 0; i < 2; i++)
	{
		/* Create a new node */
		temp = (NodePtr)malloc(sizeof(Node));
		
		/* Get the data from the file */
		fscanf(input_data, "%s", temp_id);
		fscanf(input_data, "%d", &temp_status);
		fscanf(input_data, "%lf", &temp_pay);
		fscanf(input_data, "%lf", &temp_hours);

		strcpy(temp->emp_id, temp_id);
		temp->emp_status = temp_status;
		temp->pay_rate = temp_pay;
		temp->hours_worked = temp_hours;	
		temp->link = NULL;
		
		if (temp->emp_status == 1)
		{
			temp_benefit = 
			  (0.05) * (temp->pay_rate) * (temp->hours_worked);

			temp->benefit = temp_benefit;	
		}
		else
			temp->benefit = 0;		


		/* Attach the node to the end of the list */
		if(last)
			last->link = temp;
	
		/* Update the last pointer */
		last = temp;

		/* Update the head pointer only if it is the first node */
		if (!head)
			head = temp;
	}

	/* Print the list */
	PrintList(head);
}

/*
*  Function: PrintList
*/
void
PrintList (NodePtr node_ptr)
{
	if (!node_ptr)
		printf("The list is empty.\n");

	while (node_ptr)
	{
		printf("\n%s  %d  $%.2f  %.1f  $%.2f",  	 
			node_ptr->emp_id,
			node_ptr->emp_status,
			node_ptr->pay_rate,
			node_ptr->hours_worked,
			node_ptr->benefit);

		node_ptr = node_ptr->link;
	}
	printf("\n\n");

	return;
}
