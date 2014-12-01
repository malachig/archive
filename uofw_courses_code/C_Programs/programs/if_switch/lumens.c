/* lumens.c */

/* Author: Malachi Griffith
*  Date: Oct. 25 2002
*  Purpose: To convert watts to lumens.  Use switch command.
*/

#include <stdio.h>

main()
{
	int watts;
	int lumens;
	int test;

	printf("\nEnter the watts of your bulb, enter 99 to quit > ");
	scanf("%d", &watts);

	while (watts != 99)
	{
		test = 1;

		switch (watts)
		{
		case 15:
			lumens = 125;
			break;
	
		case 25:
			lumens = 215;
			break;

		case 40:
			lumens = 500;
			break;

		case 60:
			lumens = 880;
			break;
	
		case 75:
			lumens = 1000;
			break;

		case 100:
			lumens = 1675;
			break;

		default:
			test = 0;
			
		}
		if (test == 1)		
			printf("\nThe correponding brightness is: %d\n\n",
			       lumens);
		else	
			printf("\nNot a valid entry!\n\n");
		printf("\nEnter the watts of your bulb, enter 99 to quit > ");
		scanf("%d", &watts);	
	}
}



