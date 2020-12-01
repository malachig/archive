/* Author: Malachi Griffith
*  Date:  Oct. 20 2002
*  Purpose: Finds the position of the smallest element in the subarray
*  list[first] through list[last].
*  Pre:  first < last and elements 0 through last array are defined.
*  Post: Returns the subscript k of the smallest element in the subarray;
*	i.e., list[k] <= list[i] for all i in the subarray
*/

int get_min_range(int list[], int first, int last);


/*
*  Sorts the data in array list
*  Pre: first n elements of list are defined and n >= 0
*/

void
select_sort(int list[],	/* input/output - array being sorted */
	    int n)	/* input - number of elements to sort */
{
	int fill,	/* first element in unsorted subarray */
	    temp,	/* temporary storage */
	    index_of_min; /* subscript of next smallest element */

	for (fill = 0; fill < n - 1; ++fill)
	{
		/* Find position of smallest element in unsorted subarray */
		index_of_min = get_min_range(list, fill, n-1);

		/* Exchange elements at fill and index_of_min */
		if (fill != index_of_min)
		{
			temp = list[index_of_min];
			list[index_of_min] = list[fill];
			list[fill] = temp;
		}
	}
}
