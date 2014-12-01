/* Author: Malachi Griffith
*  Date: Nov. 23 2002
*  Purpose: Illustrate the creation of an empty linked list.
*/

#include <stdio.h>

typedef struct node{
			int data;
			struct node *link;
		} Node;

typedef Node *Nodeptr;

main()
{
	Nodeptr list;
	list = NULL;
}



