/* Author: Malachi Griffith
*  Date: Dec. 11 2002  
*  Purpose: Use a queue to store the data from a datafile, consisting of
*  product_ids and the corresponding price.  This program will create the 
*  queue, print it, delete it, and search for a specified id and if it exists
*  display the price associated with it.  The data are already sorted by ID
*  and this fact will be used in the searching algorithm.
*/

#include <stdio.h>

/* Define a type for the data */
typedef struct node {
		int product_id;
		double product_price;
		struct node *link;
	} Node;

typedef Node *NodePtr;

/* Function Prototypes */
void Enqueue(NodePtr *front_ptr, NodePtr *last_ptr, int id, double price);
void PrintList(NodePtr node_ptr);
int SearchQueue(NodePtr node_ptr, int search_value, double *corr_price);
int Dequeue(NodePtr *front_ptr, NodePtr *last_ptr);
int IsEmpty(NodePtr node_ptr);
main()
{
	FILE *input_data;

	NodePtr front = NULL;
	NodePtr last = NULL;
	NodePtr temp = NULL;

	int id;
	double price;

	int search_id;
	int search_result;  /* 1 = found, 0 = not found */
	double price_found;

	int i;

	/* Create the Queue */
	input_data = fopen("stack.dat", "r");

	fscanf(input_data, "%d%lf", &id, &price);
	
	while(!feof(input_data))
	{
		Enqueue(&front, &last, id, price);
		fscanf(input_data, "%d%lf", &id, &price);
	}

	/* Print the Queue */
		PrintList(front);

	/* Search for a value specified by the user */
	printf("\nEnter the ID you would like to search for > ");
	scanf("%d", &search_id);
	
	search_result = SearchQueue(front, search_id, &price_found);

	if (search_result == 1)
		printf("\nThe Id was found, price is $%.2f\n", price_found);
	else
		printf("\nThat ID was not found, sorry dude!\n");

	/* Delete the queue */
	while(!IsEmpty(front))  /* while not returning a TRUE, ie while not
				 * empty. */
	{
		printf("The deleted node contained the ID %d\n",
			Dequeue(&front, &last));
	}

	/* Print the now empty queue */
	printf("The queue now contains:\n");
	PrintList(front);
}

/*
*  Function: Enqueue()
*/
void 
Enqueue(NodePtr *front_ptr, NodePtr *last_ptr, int id, double price)
{
	NodePtr temp;  /* points to the new node */

	/* Create a new node */
	temp = (NodePtr)malloc(sizeof(Node));

	temp->product_id = id;
	temp->product_price = price;
	temp->link = NULL;

	/* Attach the new node to the end of the list */
	if(*last_ptr)
		(*last_ptr)->link = temp;

	/* Update the last_ptr */
	(*last_ptr) = temp;

	/* Update the front pointer if it is the first node) */
	if(!(*front_ptr))
		(*front_ptr) = temp;
}

/*
*  Function: PrintList()
*/
void 
PrintList(NodePtr node_ptr)
{
	if(node_ptr == NULL)
		printf("\nThe list is empty dude!\n");

	while(node_ptr)
	{
		printf("\nProduct ID: %d, Product Price: $%.2f",
		       node_ptr->product_id,
		       node_ptr->product_price);
	
		/* Advance to next node */
		node_ptr = node_ptr->link;
	}
	printf("\n");
}

/*
*  Function: SearchQueue()
*/	
int 
SearchQueue(NodePtr node_ptr, int search_value, double *corr_price)
{
	int answer = 0;  /* 1 = found, 0 = not found */

	while ((node_ptr->product_id) <= search_value)
	{
		if(node_ptr->product_id == search_value)
		{
			answer = 1;
			*corr_price = node_ptr->product_price;
			break;
		}
	/* Advance to next node */
	node_ptr = node_ptr->link;
	}
	return(answer);	
}

/*
*  Function: Dequeue()
*/
int 
Dequeue(NodePtr *front_ptr, NodePtr *last_ptr)
{
	NodePtr temp;
	int id = 0;  /*If queue is empty this meaningless data is returned*/

	/* first check to make sure the queue is not already empty */
	if(!IsEmpty(*front_ptr))
	{
		id = (*front_ptr)->product_id;
		temp = (*front_ptr);
		(*front_ptr) = (*front_ptr)->link;
	
		if((*last_ptr) == temp)  /* If this is last node */
			(*last_ptr) = NULL;

		free(temp);
	}
	else
	 printf("\nThe queue is empty, the datum returned is meaningless.\n");
	
	return id;
}
	
/*
*  Function: IsEmpty();
*/
int 
IsEmpty(NodePtr node_ptr)
{
	int answer;

	if (node_ptr == NULL)
		answer = 1;
	else
		answer = 0;

	return(answer);
}
