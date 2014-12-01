/* Author: Malachi Griffith
*  Date:  Nov. 4 2002 
*  Purpose: A program that uses a stack to print a string in reverse.
*  Note that we have made Push and Pop more general by defining a 
*  typedef for the data type.  In this way, the same functions can be
*  used for integers and charcters, without modification.
*/

#include <stdio.h>

typedef char DataType;

typedef struct node{
	DataType data;
	struct node *link;
} Node;

typedef Node *NodePtr;

/* Function Prototypes */
int IsEmpty(NodePtr node_ptr);
void Push(NodePtr *top_ptr, DataType data);
DataType Pop(NodePtr *top_ptr);

main()
{
	NodePtr top = NULL;
	NodePtr temp;

	DataType data;

	/* Get the input line */
	printf("Please enter a line of characters:\n");

	while((data = getchar()) != '\n') 
		Push(&top, data);

	printf("\n");
	printf("The line printed in reverse is: \n");

	while(!IsEmpty(top))
	{
		data = Pop(&top);
		printf("%c", data);
	}

	printf("\n");
}

/* 
*  Function: IsEmpty();
*/
int 
IsEmpty(NodePtr node_ptr)
{
	int answer;
	answer = (node_ptr == NULL);
	return answer;
}

/*
*  Function: Push()
*/
void
Push(NodePtr *top_ptr, DataType data)
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
DataType
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
