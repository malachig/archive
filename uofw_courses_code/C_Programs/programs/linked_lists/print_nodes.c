/* Author: Malachi Griffith
*  Date: Nov. 23 2002
*  Purpose: Illustrate the basic usage of self-referential structures.
*/

#include <stdio.h>

typedef struct node{
			int data;
			struct node *link;
		   } Node;

typedef Node *NodePtr;

/* Function Prototypes */
void PrintNodes(NodePtr);

main()
{
	Node a_node;
	Node b_node;

	a_node.data = 1;
	b_node.data = 2;

	a_node.link = &b_node;
	b_node.link = NULL;

	PrintNodes(&a_node);
}

/*
*  Function: PrintNode()
*/
void
PrintNodes(NodePtr node_ptr)
{
	while(node_ptr)  /* ie. while not NULL! */
	{
		printf("%d\n", node_ptr->data);
		node_ptr = node_ptr->link;
	}
	return;
}	
