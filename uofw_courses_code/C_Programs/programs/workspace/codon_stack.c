/* Author: Malachi Griffith
*  Date: Dec. 15 2002 
*  Purpose: Create a stack to read in the data from the file codon.dat
*/

#include <stdio.h>
#include <string.h>

typedef struct node{
	char letter;
	char abbrev[5];
	char fullname[20];
	char codon[5];
	struct node *link;
		} Node;

typedef Node *NodePtr;

/* Function Prototypes */
void Push(NodePtr *top_ptr, char letter, char abbrev[], char fullname[],
	  char codon[]);
void PrintStack(NodePtr node_ptr);
void Pop(NodePtr *top_ptr);
int IsEmpty(NodePtr node_ptr);

main()
{
	NodePtr top = NULL;
	NodePtr temp = NULL;

	FILE *input_data;

	char letter;
	char abbrev[5];
	char fullname[20];
	char codon[5];

	/* Create the stack */
	input_data = fopen("codon.dat", "r");
	
	fscanf(input_data, " %c %s %s %s", &letter, abbrev, fullname, codon);

	while(!feof(input_data))
	{
		Push(&top, letter, abbrev, fullname, codon);
		fscanf(input_data, " %c %s %s %s", &letter, abbrev, 
			fullname, codon);
	}
	
	/*  Print out the stack */
	PrintStack(top);
	
	/* Delete the list */
	while(!IsEmpty(top))	
	{
		Pop(&top);
	}

	/* Print the now empty stack */
	PrintStack(top);
		
	fclose(input_data);
}

/*
*  Function: Push()
*/
void 
Push(NodePtr *top_ptr, char letter, char abbrev[], char fullname[],
     char codon[])
{
	NodePtr temp;

	/* Create a new node */
	temp = (NodePtr)malloc(sizeof(Node));

	temp->letter = letter;
	strcpy(temp->abbrev, abbrev);
	strcpy(temp->fullname, fullname);
	strcpy(temp->codon, codon);

	/* Attach the node to the list */
	temp->link = *top_ptr;

	/* Update the top pointer */
	*top_ptr = temp;
}

/*
*  Function: PrintStack()
*/ 
void 
PrintStack(NodePtr node_ptr)
{
	if(!node_ptr)
		printf("\nThe stack is empty\n");

	while(node_ptr)
	{
		printf("%c  %s  %s\t\t     %s\n", node_ptr->letter,
		       node_ptr->abbrev, node_ptr->fullname,
		       node_ptr->codon);
	
		node_ptr = node_ptr->link;
	}
}

/*
*  Function: Pop()
*/
void
Pop(NodePtr *top_ptr)
{
	NodePtr temp;
	temp = *top_ptr;
	(*top_ptr) = (*top_ptr)->link;
	free(temp);	
}

int 
IsEmpty(NodePtr node_ptr)
{
	int answer;
	answer = (node_ptr == NULL);
	return(answer);
}
