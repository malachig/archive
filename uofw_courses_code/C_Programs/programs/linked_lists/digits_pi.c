/* Author: Malachi Griffith
*  Date: Dec. 8 2002
*  Purpose: Reads the first ten thousand digits of PI from a 
*  file and inputs them into a linked list.  Each node contains 
*  one digit and a link field to the next node.
*/

#include <stdio.h>

typedef struct node{
	char digit;
	struct node *link;
		} Node;

typedef Node *NodePtr;

/* Function Prototypes */
void CreateList(NodePtr *head, NodePtr *last);
void PrintList(NodePtr head);
void DeleteList(NodePtr *head);

main()
{
	NodePtr head = NULL;
	NodePtr last = NULL;

	/* Create the list */
	CreateList(&head, &last);

	/* Display the list */
	PrintList(head);

	/* Delete the list */
	DeleteList(&head);

	/* Print the now empty list */
	printf("\nAfter freeing the memory:");
	PrintList(head);
}

/*
*  Function CreateList()
*/
void
CreateList(NodePtr *head, NodePtr *last)
{
	char ch;
	NodePtr temp;
	FILE *input_data;


	input_data = fopen("digits_pi.dat", "r");
	
	fscanf(input_data, "%c", &ch); 

	while(!feof(input_data))
	{
		/* Create the new node */
		temp = (NodePtr)malloc(sizeof(Node));
		temp->digit = ch;
		temp->link = NULL;

		/* Attach node to end of list */
		if (*last)  /* If not NULL */
			(*last)->link = temp;

		/* Update the last pointer */
		(*last) = temp;	

		/* Update the head pointer if first node */
		if (!(*head))
			(*head) = temp;
		
		fscanf(input_data, "%c", &ch); 
	}
}

/*
*  Function: PrintList()
*/
void 
PrintList(NodePtr head)
{
	if (!head)  /* or if (head == NULL) */
		printf("\nThe list is empty\n");

	while (head)
	{
		printf("%c", head->digit);
		head = head->link;
	}
}

/*
*  Function: DeleteList()
*/
void 
DeleteList(NodePtr *head)
{
	NodePtr temp;
	
	while (*head)
	{
		/* Get first node */
		temp = *head;

		/* Advance to next node before deleting */
		*head = (*head)->link;

		/* Free memory associated with node */
		free(temp);
	}
}
