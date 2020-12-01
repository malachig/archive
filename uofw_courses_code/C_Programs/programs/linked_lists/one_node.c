/* Author: Malachi Griffith
*  Date: Nov. 23 2002
*  Purpose: Creation of a linked list of one node.  Note that the printing
*  function, now renamed PrintList(), has been modified to deal with the case
*  of an empty list, where no nodes have been defined.
*/

#include <stdio.h>

typedef struct node{
			int data;
			struct node *link;
		} Node;

typedef Node *NodePtr;

/* Function Prototype */
void PrintList (NodePtr);

main()
{
	NodePtr temp;
	NodePtr list = NULL;

	temp = (NodePtr)malloc(sizeof(Node));

	temp->data = 3;
	temp->link = NULL;

	list = temp;

	PrintList(list);
}

/*
*  Function: PrintList
*/
void
PrintList(NodePtr node_ptr)
{
	/* If the node pointer is still set to NULL */
	if (!node_ptr)
		printf("The list is empty.\n");
	
	while(node_ptr)
	{
		printf("%d\n", node_ptr->data);
		node_ptr = node_ptr->link;
	}
	return;
}
