/* Author: Malachi Griffith
*  Date: Dec. 4 2002  
*  Purpose: A complete program with all neccessary functions to create
*  a stack, print it and delete it.
*/

#include <stdio.h>
#define ITER 10

typedef struct node {
	int data;
	struct node *link;
	} Node;

typedef Node *NodePtr;

/* Function Prototypes */
void Push(NodePtr *top_ptr, int data);
void PrintStack(NodePtr node_ptr);
int Pop(NodePtr *top_ptr);
int IsEmpty(NodePtr node_ptr);

main()
{
	NodePtr temp = NULL;
	NodePtr top = NULL;

	int i;
	int data;
	
	/* Create the stack */
	for(i = 0; i < ITER; i++)
		Push(&top, i + 1);

	/* Print the stack */
	PrintStack(top);

	/* Delete the stack */
	for (i = 0; i < ITER; i++)
	{
		if(!IsEmpty(top))
		{
			data = Pop(&top);
			printf("The deleted node contained %d\n", data);
		}
	}
	
	Pop(&top);

	/* Print the now-empty stack */
	PrintStack(top);
}

/*
*  Function: PrintStack()
*/
void
PrintStack(NodePtr node_ptr)
{
	if(!node_ptr)	
		printf("The stack is empty.\n");

	while(node_ptr)
	{
		printf("%d\n", node_ptr->data);
		node_ptr = node_ptr->link;
	}
	return;
}

/*
*  Function: IsEmpty()
*/
int 
IsEmpty(NodePtr node_ptr)
{
	int answer;
	answer = (node_ptr == NULL);

	return (answer);
}

/* 
*  Function: Push()
*/
void 
Push(NodePtr *top_ptr, int data)
{
	NodePtr temp;
	
	temp = (NodePtr)malloc(sizeof(Node));

	temp->data = data;

	/* Attach the new node as the top node */
	temp->link = *top_ptr;

	/* Update the top pointer */
	*top_ptr = temp;
}

/* 
*  Function: Pop()
*/
int
Pop(NodePtr *top_ptr)
{
	NodePtr temp;

	int data = 0;

	if(!IsEmpty(*top_ptr))
	{
		data = (*top_ptr)->data;
		temp = *top_ptr;
		*top_ptr = (*top_ptr)->link;
		free(temp);
	}

	else 
		printf("The stack is empty, and the datum returned is meaningless.\n");

	return (data);
}
