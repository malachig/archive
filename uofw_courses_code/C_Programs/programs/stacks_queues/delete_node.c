/* Author: Malachi Griffith
*  Date: Nov. 3 2002 
*  Purpose: This program deletes the top node of the stack, and prints the
*  rest of the stack.
*/

#include <stdio.h>
#define ITER 10

typedef struct node{
	int data;
	struct node *link;
	} Node;

typedef Node *NodePtr;

/* Function Prototypes */
void PrintNodes(NodePtr);
void Push(NodePtr *top_ptr, int data);

main()
{
	NodePtr top = NULL;
	NodePtr temp;
	int i;

	/* create the stack */
	for(i = 0; i < ITER; i++)
		Push(&top, i + 1);

	/* delete the top node */
	if (top)
	{
		printf("Top node datum is %d,\n", top->data);
		
		/* get top node */	
		temp = top;
	
		/* advance to next node before deleting */
		top = top->link;
	
		/* delete top node by freeing temp */
		free(temp);

		printf("This node has been deleted.\n");
	}
	else 
		printf("The stack is empty.\n");

	/* print the stack */
	PrintNodes(top);
}

/*
*  Function: PrintNodes();
*/
void
PrintNodes(NodePtr node_ptr)
{
	while(node_ptr)
	{
		printf("%d\n", node_ptr->data);
		node_ptr = node_ptr->link;
	}
}

void
Push(NodePtr *top_ptr, int data)
{
	NodePtr temp;

	temp = (NodePtr)malloc(sizeof(Node));
	temp->data = data;

	/* attach the node as the top node */
	temp->link = *top_ptr;

	/* update the top pointer */
	*top_ptr = temp;
}
