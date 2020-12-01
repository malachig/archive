/* Author: Malachi Griffith
*  Date: Nov. 24 2002 
*  Purpose: This programs uses a linked list to diplay data from 
*  class 1 sales only.  No sorting is required.
*/

#include <stdio.h>

#define COM_RATE1 0.045  /* Rate of commission for class 1 property sales */

/* Define a structure to store property sale records.  Include a link
*  field to allow use of self-referential structures. */
typedef struct node {
			int property_class;  	  /* 1, 2 or 3 */
			double selling_price;     /* Sale price */
			double commission_earned; /* Price x Rate */
			struct node *link;	  /* to create linked list */
	} Node;  /* Name of this structure */

/* Now create Pointer variable for the structure 'Node' */
typedef Node *NodePtr; 

/* Function Prototypes */
void make_list(NodePtr *head, NodePtr *last);
void display_list(NodePtr node_ptr);
void delete_list(NodePtr *head);

main()
{
	NodePtr head = NULL;  /* Initialize head to NULL */
	NodePtr last = NULL;  /* Initialize last to NULL */

	/* Make the linked list - User Defined Funtion */
	make_list(&head, &last);

	/* Display the list - User Defined Function */
	display_list(head);

	/* Delete the List - free memory space - User Defined Function */
	delete_list(&head);
	
	/* Verify the list has been successfully deleted */
	printf("After deletion of linked list:\n");
	display_list(head);
}

/*
*  Function: make_list()
*  Creates a linked list of nodes where each node is a structure containing
*  a sales record.  Only entries from the file for class 1 property sales
*  are added to the linked list!.
*  Pre: The function receives no data from main.
*  Post: 'head' and 'last' in main are assigned values. 
*/
void make_list(NodePtr *head, NodePtr *last)
{
	FILE *input_data;  /* File pointer variable */

	NodePtr temp;   /* Points to each new node */
	int class; 	/* Temp local variable for the property class */
	double price;   /* Temp local variable for the sale price */
	double commission;  /* Temp local variable for commission earned */

	/* Open the input file for reading */
	input_data = fopen("assign4.dat", "r");

	/* Get a priming set of data from the input file */	
	fscanf(input_data, "%d%lf", &class, &price);

	/* Continue adding nodes until EOF as long as record is class 1 */
	while(!feof(input_data))
	{
		/* Create the new node - using dynamic memory allocation */
		temp = (NodePtr)malloc(sizeof(Node));
		
		/* Attach the node if the class is 1 */	
		if (class == 1)
		{
		   /*  Set each value in the temp structure */
		   temp->property_class = class;
		   temp->selling_price = price;
		   
		   commission = (temp->selling_price) * COM_RATE1;
		  
		   temp->commission_earned = commission;
		   temp->link = NULL; /* Initialize link to NULL */

		   /* Attach the node to the end of the list */
		   if (*last) 
			(*last)->link = temp;

		   /* Update the last pointer */
		   (*last) = temp;

		   /* Update the head pointer only if it is the first node */
		   if (!(*head))
			(*head) = temp;
		}
		/* Get the next set of data for analysis */
		fscanf(input_data, "%d%lf", &class, &price);
	}
	fclose(input_data);
}

/*
*  Function: display_list()
*  Simply receives the pointer to the head node and traverses the linked
*  list, displaying the sales records to the screen as it goes.
*  Pre: node_ptr is a pointer variable of NodePtr type.  
*  Post: Displays records to screen.
*/
void 
display_list(NodePtr node_ptr)
{
	double total_commission = 0;  /* Initialize accumulator to 0 */

	/* Check to make sure the list is not empty, i.e. not NULL! */
	if (!node_ptr)
	{
		printf("The list is empty.\n");
		/* If it is empty, exit gracefully */
		exit(1);
	}
	
	/* Display a heading */
	printf("\nProperty Class\t\tSelling Price\t\tCommission Earned\n");

	while (node_ptr)
	{
		printf("     %d\t\t        $%.2f\t\t    $%.2f\n",
			node_ptr->property_class,
			node_ptr->selling_price,
			node_ptr->commission_earned);
		
		/* Tally the commission earned for this class */	
		total_commission += node_ptr->commission_earned;

		/* Advance to next node in list */
		node_ptr = node_ptr->link;
	}
	printf("\nTotal commission earned was $%.2f\n\n", total_commission);

	return;
}

/*
*  Function: delete_list()
*  Simply frees the memory used by the linked list after it is no longer 
*  needed.
*  Pre:  The location of first node is specified so that the function can
*  traverse the linked list and delete each node as it goes.
*  Post: Updates head until head reaches the last link which is NULL.
*/
void
delete_list(NodePtr *head)
{
	NodePtr temp;

	while (*head)
	{
	   temp = *head;
	   *head = (*head)->link; /*Advance to next node before freeing*/
	   free(temp);
	}
}
