/* noise_db.c */

/* Author: Malachi Griffith
*  Date:  Oct. 2 2002
*  Purpose: To illustrate the use of if-else statements for 
*	    multiple alternative decisions.
*/

#include <stdio.h>

main()
{
	/* Display perception of noise loudness */

	int noise_db;

	printf("Enter the noise level in decibels (0-150) > ");
	scanf(" %d", &noise_db);
	if (noise_db <= 50)
		printf("%d-decibel noise is quiet.\n", noise_db);
	else if (noise_db <= 70)
		printf("%d-decibel noise is intrusive.\n", noise_db);
	else if (noise_db <=90)
		printf("%d-decibel noise is annoying.\n", noise_db);
	else if (noise_db <= 110)
		printf("%d-decibel noise is very annoying.\n", noise_db);
	else
		printf("%d-decibel noise is uncomfortable.\n", noise_db);
}



