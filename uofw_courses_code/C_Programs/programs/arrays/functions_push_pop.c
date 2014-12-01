/* Author: Malachi Griffith
*  Date: Oct. 20 2002
*  Purpose: To illustrate basic stacking functions.
*/

/* Function: push */

void
push(char stack[],  /* input/output - the stack */
     char item,     /* input - data being pushed onto the stack */
     int *top,	    /* input/output - pointer to top of stack */
     int max_size)  /* input - maximum size of stack */

{
	if (*top < max_size - 1)
	{
	++(*top);
	stack[*top] = item;
	}
}


/* Function: pop */
char
pop(char stack[],	/* input/output - the stack */
    int *top)		/* input/output - pointer to top of stack */
{
	char item;	/* value popped of the stack */
		
	if (*top >= o)
	{
		item = stack[*top];
		--(*top);
	}
	else
	{
		item = STACK_EMPTY;
	}
	
	return(item);
}

