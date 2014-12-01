/* Author: Malachi Griffith
*  Date: Oct. 19 2002
*  Purpose: Returns the largest of the first n values in array list
*  Pre:  First n elements of array list are defined and n > 0.
*/

int 
get_max(const int list[],  /* input - list of n integers */
	int n)		   /* input - number of list elements to examine */

/* NOTE: the reserved term 'const' specifies that the array variable 
*  declared is is striclty an input parameter and will not be modified 
*  by the function */

{
	int i,
	cur_large;	/* Largest value so far */

	/* initial array element is largest so far. */
	cur_large = list[0];
	
	/* Compare each remaining list element to the largest so far; 
	*  save the larger */

	for (i = 1; i < n; ++i)
		if (list[i] > curr_large)
			cur_large = list[i];

	return (cur_large);
}

/* NOTE: To call this function you would use something like,
*		
*	x_large = get_max(x, 5);
*
*	Where x is an array which is searched for its largest value up to 
*	element 5 in the array.
*/

