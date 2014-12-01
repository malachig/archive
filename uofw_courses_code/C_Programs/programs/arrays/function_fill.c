/* Author: Malachi Griffith
*  Date: Oct. 19 2002
*  Purpose: A function that sets all elements of its array parameter
*  to in_value.
*  Pre: n and in_value are defined.
*  Post: list[i] = in_value, for 0 <= i < n.
*/

void
fill_array (int list[],		/* output - list of n integers */
	    int n,		/* input - number of list elements */
	    int in_value)	/* input - initial value */

{
	int i;		/* array subscript and loop control */

	for (i = 0; i < n; ++i)
		list[i] = in_value;
}

/* NOTE:  To call function fill_array, you must specify the actual array
*  argument, the number of array elements, and the value to be stored
*  in the array.  If y is an array with ten type int elements, the function
*  call,
*  		fill_array(y, 10, num);
*
*  stores the value of num in the ten elements of array y.
*/


