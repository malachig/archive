/* Author: Malachi Griffith
*  Date: Oct. 20 2002
*  Purpose: Function that searches for target item in first element of array
*  arr and Returns index of target or NOT_FOUND.
*  Pre: target and first n elements of array arr are defined n>=0.
*/

#define NOT_FOUND -1

int
search (const int arr[],	/* input- array to search */
	int 	target,		/* input value searched for */
	int 	n 		/* input - number of elements to search */
{
	int i,
	found = 0,	/* whether or not the target has been found */
	where;		/* index where target found or NOT_FOUND */

	/* Compares each element to target */
	i = 0;
	while (!found && i < n)
	{
		if (arr[i] == target)
			found = 1;
		else
			++i;
	}

	/* Returns index of element matching target or NOT_FOUND */
	if (found)
		where = i;
	else 
		where = NOT_FOUND;
	
	return(where);
}
