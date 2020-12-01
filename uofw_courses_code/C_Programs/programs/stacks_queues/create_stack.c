/* Author: Malachi Griffith
*  Date: Dec. 4 2002  
*  Purpose: This program creates a stack one node at a time, and then
*  prints it.
*/

#include <stdio.h>
#define ITER 10

typedef struct node{
	int data;
	struct node *link;
	} Node;

typedef Node *NodePtr;

void PrintNodes(NodePtr);

main()
{
	NodePtr temp;   /* used to point at new nodes */
	NodePtr top = NULL;  /* used to point at the first node */

	int i;

	for(i = 0; i < ITER; i++)
	{
		/* Create a new node */
		temp = (NodePtr)malloc(sizeof(Node));
		temp->data = i + 1;	
		
		/* Attach the new node as top node */
		temp->link = top;

		/* update the top pointer */
		top = temp;
	}

	/* Print the stack */
	PrintNodes(top);
}

/*
*  Function: PrintNodes
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


