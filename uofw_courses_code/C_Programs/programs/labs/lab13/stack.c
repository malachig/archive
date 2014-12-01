/* Author: Malachi Griffith
*  Date: Nov. 27 2002 
*  Purpose: Create a stack to contain the product information in the 
*  data file "stack.dat".  The program will print the stack and delete the 
*  stack at the end of the program.  It will also search for a product of
*  a particular price and print out its' ID.
*/

#include <stdio.h>

typedef struct node{
			int product_id;
			double product_price;
			struct node *link;
		} Node;

typedef Node *NodePtr;

void Push(NodePtr *top_ptr, int id, double price);
void PrintStack(NodePtr node_ptr);
int Pop(NodePtr *top_ptr);
int IsEmpty(NodePtr node_ptr);
int SearchStack(NodePtr node_ptr, double search_value);

main()
{
	FILE *input_data;

	NodePtr temp = NULL;
	NodePtr top = NULL;

	int i;
	int id;
	double price;
	int size = 0;
	double search_value;
	int search_id;

	input_data = fopen("stack.dat", "r");

	fscanf(input_data, "%d%lf", &id, &price);
	/* Create the stack */
	while (!feof(input_data))
	{	
		Push(&top, id, price);
		size++;
		fscanf(input_data, "%d%lf", &id, &price);
	
	}

	/* Print the stack */
	PrintStack(top);
	
	/* Search the stack for an ID number by price */
	printf("\nEnter the price you wish to search for > ");
	scanf("%lf", &search_value);

	search_id = SearchStack(top, search_value);
	
	if (search_id != 0)
	  printf("\nThe Product ID for that price is %d\n", search_id);
	else
	  printf("\nThere is no product with that price\n");

	/* Delete the stack */
	for(i = 0; i < size; i++)
	{
		if(!IsEmpty(top))
		{
			id = Pop(&top);
			printf("The deleted node contained product %d\n", id);
		}
	}

	/* Print the now-empty stack */
	PrintStack(top);
}

/*
*  Function: PrintStack
*/
void
PrintStack(NodePtr node_ptr)
{
	if (!node_ptr)
		printf("The stack is empty.\n");

	while(node_ptr)
	{
		printf("%d\t$%.2f\n", node_ptr->product_id, 
			node_ptr->product_price);
		
		node_ptr = node_ptr->link;
	}
	return;
}

/*
*  Function: IsEmpty
*/
int
IsEmpty(NodePtr node_ptr)
{
	int answer;
	
	answer = (node_ptr == NULL);

	return (answer);
}

/*
*  Function Push()
*/
void
Push(NodePtr *top_ptr, int id, double price)
{
	NodePtr temp;

	temp = (NodePtr)malloc(sizeof(Node));

	temp->product_id = id;
	temp->product_price = price;

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

	int id = 0;

	if(!IsEmpty(*top_ptr))
	{
		id = (*top_ptr)->product_id;

		temp = *top_ptr;
		*top_ptr = (*top_ptr)->link;
		free(temp);
	}

	else
		printf("The stack is empty, and the datum returned is meaningless.\n");

	return (id);
}
/*
*  Function: SearchStack()
*/
int 
SearchStack(NodePtr node_ptr, double search_value)
{
	int id_to_return = 0;
	double temp;
	
	while(node_ptr)
	{
		temp = node_ptr->product_price;
		if (temp == search_value)
		{
			id_to_return = node_ptr->product_id;
			break;
		}

		node_ptr = node_ptr->link;
	}
	return(id_to_return);
}


