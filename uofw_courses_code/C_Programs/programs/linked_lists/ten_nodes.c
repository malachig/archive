/* Author: Malachi Griffith
*  Date:  Nov. 23 2002 
*  Purpose: Creates a linked list of ten nodes, each of type Node.
*  The first node will contain the integer value 1, the second 2, etc.,
*  with the last node containing the integer value 10.
*/

#include <stdio.h>
#define ITER 10


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
	NodePtr head = NULL;
	NodePtr last = NULL;

	int i;
	int key;

	for (i = 0; i < ITER; i++)
	{
		/* Create the new Node */
		temp = (NodePtr)malloc(sizeof(Node));

		temp->data = i + 1;
		temp->link = NULL;
	
		/* Attach the node to the end of the list */
		
		if (last)
			last->link = temp;

		/* Update the last pointer */
		
		last = temp;

		/* Update the head pointer only if it is the first node */

		if (!head)
			head = temp;
	}
	
	/* Print the list */
	PrintList (head);
}

/*
* Function: PrintList()
*/
void
PrintList(NodePtr node_ptr)
{
	if (!node_ptr)
		printf("The list is empty.\n");

	while (node_ptr)
	{
		printf("%d\n", node_ptr->data);
		
		node_ptr = node_ptr->link;
	}
	return;
}

/*  DISCUSSION AND EXPLANATION 
*  
*  1.  Let 'head' be a pointer that always points at the first node, analogous
*      to an array name pointer.
*
*  2.  Let 'temp' be a pointer that always points to a newly created node:
*	
*	a) temp = (NodePtr)malloc(sizeof(Node));
*	b) temp -> data = i + 1;
*	c) temp -> link = NULL;
*
*  3.  Let 'last' be a pointer that always points to the last node of the 
*      list, i.e., 'last' has a value different from NULL if there is at
*      least one node in the list.  We attach the new node to the end of the 
*      list by making the pointer field of the current last node point to 
*      the address of the new node:
*	
*	if (last)
*		last -> link = temp;
*	
*  4.  We update the 'last' pointer to point to the node we have just added,
*      which is now the last node in the list:
*
*	last = temp;
*
*  5.  If the newly created node is the first node, i.e., if 'head' is still
*      at its initial value of NULL, then we set 'head' to point to the new
*      node:
*
*	if (!head)
*		head = temp;
*
*  We can repeat steps 2-5 as many times as we wich to create nodes.
*/



