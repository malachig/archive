/* Author: Malachi Griffith
*  Date: Dec. 8 2002 
*  Purpose: Simply a function to display the list
*  pointed to by headp, a pointer to the first node in 
*  a list.
*/

typedef struct list_node_s {
	int digit;
	struct list_node_s *restp;
		} list_node_t;

/* USING RECURSION */
void
print_list(list_node_t *headp)
{
	if (headp == NULL)
		printf("\n");
	else
	{
		printf("%d", headp->digit);
		print_list(headp->restp);
	}
}

/* Without Recursion */
void
print_list(list_node_t *headp)
{
	if(head == NULL)
		printf("\n");

	while(head != NULL)
	{
		printf("%d", headp->digit);
		headp = headp->restp;
	}
