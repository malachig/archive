/* Author: Malachi Griffith
*  Date: Dec. 8 2002  
*  Purpose: Forms a linked list of an input list of integers terminated
*  by SENT
*/

list_node_t *
get_list(void)
{
	int data;
	list_node_t *ansp,
		    *to_fillp,  /* pointer to last node in list whose
				   restp component is unfilled */
		    *newp;	/* pointer to newly allocated node */

	/* Builds first node, if there is one */
	scanf("%d", &data);
	if (data == SENT)
		ansp = NULL;
	else
	{
		ansp = (list_node_t *)malloc(sizeof(list_node_t));
		ansp->digit = data;
		to_fillp = ansp;

		/* Continues building list be creating a node on each
		iteration and storing its pointer in the restp component
		of the node accessed through to fillp */

		for (scanf("%d", &data);
		     data != SENT;
		     scanf("%d", &data))
		{
			newp = (list_node_t *)malloc(sizeof(list_node_t));
			newp->digit = data;
			to_fillp->restp = newp;
			to_fillp = newp;
		}

		/* Stores NULL in final node's restp component */
		to_fillp->restp = NULL;
	}
	return(ansp);
}

	








