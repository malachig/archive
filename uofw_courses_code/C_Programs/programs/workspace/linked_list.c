#include <stdio.h>

typedef struct node{
	int number;
	struct node *link;
	} Node;

typedef Node *NodePtr;

/* Function Prototypes */
void CreateList(NodePtr *head, NodePtr *last);
void PrintList(NodePtr node_ptr);
int DeleteList(NodePtr *head);

main()
{
	NodePtr head = NULL;
	NodePtr last = NULL;
	int last_data;

	/* Create the list */
	CreateList(&head, &last);

	/* Print the list */
	PrintList(head);

	/* Delete the list */
	last_data = DeleteList(&head);
	
	/* Print the empty list */
	PrintList(head);
	printf("\nThe last data point to be deleted was %d\n", last_data);
}

void
CreateList(NodePtr *head, NodePtr *last)
{
	FILE *input_data;
	NodePtr temp;
	int data;

	input_data = fopen("numbers.dat", "r");
	
	temp = (NodePtr)malloc(sizeof(Node));
	fscanf(input_data, "%d", &data);

	while(!feof(input_data))
	{
		/* Create the new node */
		temp->number = data;
		temp->link = NULL;

		/* Attach the new node */
		if (*last)	
			(*last)->link = temp;

		/* Update the last pointer */
		(*last) = temp;

		/* Update the head pointer if first node */
		if (!(*head))
			(*head) = temp;

		temp = (NodePtr)malloc(sizeof(Node));
		fscanf(input_data, "%d", &data);
	}
}

void 
PrintList(NodePtr node_ptr)
{
	if (!node_ptr)
		printf("\nThe list is empty\n");

	while(node_ptr)
	{
		printf("\n%d", node_ptr->number);
		node_ptr = node_ptr->link;
	}
}

int 
DeleteList(NodePtr *head)
{
	NodePtr temp;
	int data;

	while(*head)
	{
		temp = (*head);
		data = temp->number;
		(*head) = (*head)->link;
		free(temp);
	}
	return(data);
}
