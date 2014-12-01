/* switch.c */

/* Author: Malachi Griffith
*  Date: Oct. 2 2002
*  Purpose:  Illustrate the basic usage of the switch command
*/

#include <stdio.h>

main()
{
	char class;
	printf("Enter the type of ship, (B, C, D, or F) > ");
	scanf(" %c", &class);

	switch (class)
	{
	case 'B':
	case 'b':
		printf("\nBattleship\n");
		break;

	case 'C':
	case 'c':
		printf("\nCruiser\n");
		break;

	case 'D':
	case 'd':
		printf("\nBattleship\n");
		break;
	
	case 'F':
	case 'f':
		printf("\nFrigate\n");
		break;
	
	default:
		printf("\nUnknown ship class %c\n");	
	}
}



