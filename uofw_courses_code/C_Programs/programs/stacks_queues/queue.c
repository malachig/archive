/* Author: Malachi Griffith
*  Date: Dec. 7 2002 
*  Purpose: A queue is a FIFO (first in first out) data structure.
*/

#include <stdio.h>
#define ITER 3

/* To make it easier to change data type at a later time */
typedef int DataType;

typedef struct node{
	DataType data;
	struct node *link;
	} Node;

typedef Node *NodePtr;

/* Function Prototypes */
void Enqueue(NodePtr *front_ptr, NodePtr *last_ptr, DataType data);
void PrintList(NodePtr node_ptr);
DataType Dequeue(NodePtr *front_ptr, NodePtr *last_ptr);
int IsEmpty(NodePtr node_ptr);

main(void)
{
	NodePtr temp = NULL;
	NodePtr front = NULL;
	NodePtr last = NULL;

	DataType data;

	int i;

	/* Create the queue */
	for (i = 0; i < ITER; i++)
		Enqueue(&front, &last, i + 1);

	/* Print the queue */
	PrintList(front);

	/* Delete a node */
	if(!IsEmpty(front))
		printf("The deleted node contained %d\n",
			Dequeue(&front, &last));

	/* Print the queue */
	printf("The queue now contains:\n");
	PrintList(front);
}

/*
*  Function: PrintList()
*/
void 
PrintList(NodePtr node_ptr)
{
	if(!node_ptr)
		printf("The list is empty.\n");

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

	return answer;
}

/*
*  Function: Enqueue()
*/
void
Enqueue(NodePtr *front_ptr, NodePtr *last_ptr, DataType data)
{
	NodePtr temp;  /* points to the new node */

	/* Create a new node */
	temp = (NodePtr)malloc(sizeof(Node));

	temp->data = data;
	temp->link = NULL;

	/* Attach the new node to the end of the list */
	if (*last_ptr)	/* If last is NOT NULL */
		(*last_ptr)->link = temp;

	/* Update the last pointer */
	*last_ptr = temp;

	/* Update the front pointer if it is the first node */
	if(!(*front_ptr))  /* If front IS NULL */
		(*front_ptr) = temp;
}

/*
* Function: Dequeue()
*/
DataType
Dequeue(NodePtr *front_ptr, NodePtr *last_ptr)
{
	NodePtr temp;
	DataType data = 0;

	if (!IsEmpty(*front_ptr))
	{
		data = (*front_ptr)->data;
		temp = *front_ptr;
		*front_ptr = (*front_ptr)->link;
	
		if((*last_ptr) == temp)  /* last node? */
			(*last_ptr) = NULL;
	
		free(temp);
	}
	else
		printf("The queue is empty, and the datum returned is 
			meaningless.\n");

	return data;
}
 





