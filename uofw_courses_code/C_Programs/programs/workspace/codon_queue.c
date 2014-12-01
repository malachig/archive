/* Author: Malachi Griffith
*  Date: Dec. 15 2002 
*  Purpose: Use a Queue to store and display codon data.
*/

#include <stdio.h>
#include <string.h>

typedef struct node{
		char letter;
		char abbrev[5];
		char fullname[20];
		char codon[5];
		struct node *link;
	} Node;

typedef Node *NodePtr;

/* Function Prototypes */
void Enqueue(NodePtr *front_ptr, NodePtr *last_ptr, char letter, 
	     char abbrev[], char fullname[], char codon[]);
void PrintQueue(NodePtr node_ptr);
void Dequeue(NodePtr *front_ptr, NodePtr *last_ptr);
int IsEmpty(NodePtr node_ptr);

main()
{
	NodePtr temp = NULL;
	NodePtr front = NULL;
	NodePtr last = NULL;
	
	char letter;
	char abbrev[5];
	char fullname[20];
	char codon[5];
	
	FILE *input_data;
	
	input_data = fopen("codon.dat", "r");

	/* Create the list */
	fscanf(input_data, " %c %s %s %s", &letter, abbrev, fullname, codon);

	while(!feof(input_data))
	{
		Enqueue(&front, &last, letter, abbrev, fullname, codon);
		fscanf(input_data, " %c %s %s %s", &letter, abbrev, 
			fullname, codon);
	}

	/* Print the list */
	PrintQueue(front);
	
	/* Delete the queue */
	while (!IsEmpty)
	{
		Dequeue(&front, &last);
	}
}

/*
*  Function: Enqueue()
*/
void 
Enqueue(NodePtr *front_ptr, NodePtr *last_ptr, char letter, 
	     char abbrev[], char fullname[], char codon[])
{
	NodePtr temp;

	/* Create the new node */
	temp = (NodePtr)malloc(sizeof(Node));
	temp->letter = letter;
	strcpy(temp->abbrev, abbrev);
	strcpy(temp->fullname, fullname);
	strcpy(temp->codon, codon);
	temp->link = NULL;

	/* Attach the new node */
	if (*last_ptr)
		(*last_ptr)->link = temp;

	/* Update the last pointer */
	(*last_ptr) = temp;

	/* Update the first pointer if it is the first node */
	if (!(*front_ptr))
		(*front_ptr) = temp;
}

/*
*  Function: PrintQueue()
*/
void 
PrintQueue(NodePtr node_ptr)
{
	if (!node_ptr)
		printf("\nThe queue is empty!\n");

	while (node_ptr)
	{
		printf("%c  %s  %s  %s\n", node_ptr->letter,
		      node_ptr->abbrev, node_ptr->fullname, node_ptr->codon);
		node_ptr = node_ptr->link;
	}
}

/*
*  Function: Dequeue()
*/
void 
Dequeue(NodePtr *front_ptr, NodePtr *last_ptr)
{
	NodePtr temp;
	
	temp = (*front_ptr);
	(*front_ptr) = (*front_ptr)->link;

	if ((*last_ptr == temp))
		(*last_ptr) = NULL;

	free(temp);
}

/*
*  Function: IsEmpty()
*/
int 
IsEmpty(NodePtr node_ptr)
{
	int answer;
	answer = (node_ptr == NULL);
	
	return(answer);
}


