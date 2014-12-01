/* Author: Malachi Griffith
*  Date: Dec. 15 2002 
*  Purpose: Use a linked list to store and search through codon and 
*  amino acid data from the file codon.dat
*/

#include <stdio.h>
#include <string.h>
#define MAX_NAME 20

typedef struct node{
	char letter;
	char abbrev[5];
	char fullname[MAX_NAME];
	char codon[5];
	struct node *link;
		} Node; 

typedef Node *NodePtr;

/* Function Prototypes */
void CreateList(NodePtr *head, NodePtr *last);
void PrintList(NodePtr node_ptr);
NodePtr SearchList(NodePtr node_ptr, char codon[]);
void DeleteList(NodePtr *head);

main()
{
	NodePtr head = NULL;
	NodePtr last = NULL;
	
	char search_codon[5];
	NodePtr answer;

	/* Create the List */
	CreateList(&head, &last);

	/* Display the list */
	PrintList(head);

	/* Search the list */
	printf("\nEnter the codon (uppercase) to search for > ");
	scanf("%s", search_codon);

	answer = SearchList(head, search_codon);
	printf("\nThe corrsponding amino acid is %c, %s, %s\n", 
		answer->letter, answer->abbrev, answer->fullname);

	/* Delete the list */
	DeleteList(&head);

	/* Print the now empty list */
	PrintList(head);
}

/* 
*  Function: CreateList()
*/
void
CreateList(NodePtr *head, NodePtr *last)
{
	NodePtr temp;
	FILE *input_data;

	input_data = fopen("codon.dat", "r");
	
	/* Create the new node */
	temp = (NodePtr)malloc(sizeof(Node));
	fscanf(input_data, " %c %s %s %s", &temp->letter, temp->abbrev,
	       temp->fullname, temp->codon);
	temp->link = NULL;

	while(!feof(input_data))
	{
	  /* Attach the new node */
	  if (*last)
	    (*last)->link = temp;

	  /* Update the last pointer */
	  *last = temp;

	  /* Update the head pointer if it is the first */
	  if (!(*head))
	    (*head) = temp;

	  /* Create the new node */
	  temp = (NodePtr)malloc(sizeof(Node));
	  fscanf(input_data, " %c %s %s %s", &temp->letter, temp->abbrev,
	         temp->fullname, temp->codon);
	  temp->link = NULL;
	}
}

/*
*  Function PrintList()
*/
void
PrintList(NodePtr node_ptr)
{
	printf("\n");

	if (!(node_ptr))
		printf("\nThe list is empty\n");

	while(node_ptr)
	{
		printf("%c  %s  %s   \t\t%s\n", node_ptr->letter,
			node_ptr->abbrev, node_ptr->fullname,
			node_ptr->codon);
		node_ptr = node_ptr->link;
	}
}

/*
*  Function: SearchList()
*/	
NodePtr 
SearchList(NodePtr node_ptr, char codon[])
{
	int value;
	NodePtr temp;

	while(node_ptr)
	{
		value = strcmp(codon, node_ptr->codon);
		if (value == 0)
		{
			temp = node_ptr;
			break;
		}
		node_ptr = node_ptr->link;
	}
	return(temp);	
}

/*
*  Function: DeleteList
*/
void
DeleteList(NodePtr *head)
{
	NodePtr temp;

	while((*head))
	{
		temp = (*head);
		(*head) = (*head)->link;
		free(temp);
	}
}



