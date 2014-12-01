/* Author: Malachi Griffith
*  Date: Nov. 23 2002
*  Purpose: Search a linked list for a key value, using a search function.
*  This function will return NULL if the key is not found, otherwise it 
*  will return a pointer to the node containing the key.  Note: care must 
*  be taken never to follow a NULL pointer.
*/

#include <stdio.h>
#define ITER 10

typedef struct node{
			int data;
			struct node *link;
		} Node;

typedef Node *NodePtr;

/* Function Prototypes */
void PrintList(NodePtr);
NodePtr SearchList(NodePtr, int);

main()
{
	NodePtr temp;
	NodePtr head = NULL;
	NodePtr last = NULL;

	int i;
	int key;

	for (i = 0; i < ITER; i++)
	{
		/* Create the new node */
		
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
	PrintList(head);

	/* Search the list */
	while(!feof(stdin))
	{
		scanf("%d\n", &key);

		temp = SearchList(head, key);

		if (temp)
			printf("Key %d found.\n", key);
		else
			printf("Key %d not found.\n", key);
	}
	
	/* Delete the list to free memory used.  To ensure that we dont lose
	*  lost the ability to access a node before we have a chance to free
	*  it we assign a temporary pointer to point to the memory to be 
	*  freed, and move 'head' on to the next node before calling free()*/

	while (head)
	{
		temp = head;
		head = head->link;
		free(temp);
	}
	
	return;
}

/*
*  Function: PrintList()
*/
void
PrintList(NodePtr node_ptr)
{
	if(!node_ptr)
		printf("The list is empty.\n");

	while (node_ptr)
	{
		printf("%d\n", node_ptr->data);
		node_ptr = node_ptr->link;
	}
	return;
}

/* 
*  Function: SearchList()
*/
NodePtr
SearchList(NodePtr list, int key)
{
	NodePtr temp = NULL;

	while (list)
	{
		if (list->data == key)
		{
			temp = list;
			break;
		}
		else
			list = list->link;
	}
	return(temp);
}
