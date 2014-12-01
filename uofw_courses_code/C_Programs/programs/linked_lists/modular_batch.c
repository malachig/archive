/* Author: Malachi Griffith
*  Date: Nov. 23 2002 
*  Purpose: A modular program which creates a linked list, displays 
*  the list of data, and searches it for key values entered by the user.
*/

#include <stdio.h>
#define ITER 10

typedef struct node{
			int data;
			struct node *link;
		} Node;

typedef Node *NodePtr;

/* Function Prototypes */
void CreateList(NodePtr *, NodePtr *);
void DeleteList(NodePtr *);
void PrintList(NodePtr);
NodePtr SearchList(NodePtr, int);
FILE * SafeFopen(char *, char *);

main()
{
	NodePtr temp = NULL;
	NodePtr head = NULL;
	NodePtr last = NULL;

	int key;

	/* Create the list */
	CreateList(&head, &last);

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
	/* Delete the list */
	DeleteList(&head);

	/* Now Print the empty list */
	printf("\nAfter Freeing Memory used by List:\n");
	PrintList(head);
	
	printf("\n\n");
}

/*
*  Function: PrintList()
*/
void
PrintList(NodePtr node_ptr)
{
	if (!node_ptr)
		printf("The list is empty.\n");

	while (node_ptr)
	{
		printf("%d\n", node_ptr->data);
		node_ptr = node_ptr->link;  /* Advance to next node */
	}
	return;
}

/*
* Function: SearchList
*/
NodePtr
SearchList(NodePtr list, int key)
{
	NodePtr temp = NULL;

	while(list)
	{
		if (list -> data == key)
		{
			temp = list;
			break;
		}
		else
			list = list -> link;
	}
	return temp;
}

/* 
*  Function: CreateList()
*/
void
CreateList(NodePtr *head, NodePtr *last)
{
	FILE *input_file;

	NodePtr temp;	/* points to new node */

	int dum;	/* Reads in data from a file */

	input_file = SafeFopen("t1_input", "r");

	while(!feof(input_file))
	{
		/* Create the new node */
		temp = (NodePtr)malloc(sizeof(Node));
	
		fscanf(input_file, "%d\n", &dum);

		temp->data = dum;
		temp->link = NULL;

		/* Attach the node to the end of the list */
		if (*last)
			(*last)->link = temp;

		/* Update the last pointer */
			(*last) = temp;

		/* Update the head pointer only if it is the first node */
			if (!(*head))
				(*head) = temp;
	}
	fclose(input_file);
}

/*
* Function: SafeFopen
*/
FILE * SafeFopen(char *file_name, char *mode)
{
	FILE *file_ptr;

	if ((file_ptr = fopen (file_name, mode)) == NULL)
	{
		fprintf(stderr, "Cannot open file %s.\n", file_name);
		exit(1);
	}

	return(file_ptr);
}

/* 
*  Function: DeleteList()
*/
void
DeleteList(NodePtr *head)
{
	NodePtr temp;

	while (*head)
	{
		temp = *head;
		*head = (*head)->link;
		free(temp);
	}
}

/*  DISCUSSION
*  
*   A.)  In CreateList(), the addresses of the pointers 'head' and 'last'
*   	 are passed to the function, because we need to change their content
*     	 (the addresses which they point to) outside the function.
*
*   B.)	 In SearchList() and PrintList(), we only pass the pointer 'head',
*   	 in a manner analagous to passing the name of an array to a function.
*
*   C.)  In DeleteList(), only the address of head is passed.
*/

		





