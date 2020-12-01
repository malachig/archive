/*log_operators.c*/
/*
*  Author: Malachi Griffith
*  Date: Sept. 17 2002
*  Purpose: Ilustrate the use of logical operators
*	"&&" is like saying AND - Conjunction operator.
*	"||" is like saying OR - Disjunction operator.
*	"!" is like saying NOT - Negation operator.

/* A program to illustrate conjunction and disjunction operators */

#include <stdio.h>

main()
{
	int i = 10;
	int j = 0;
	int k = 15;

	printf("Conjoining i and j gives: %d\n", i && j);
	printf("Disjoining i and j gives: %d\n", i || j);
	printf("Conjoining i and k gives: %d\n", i && k);
	printf("Disjoining i and k gives: %d\n", i || k);

/* 
*  Note: 0 = FALSE and any NON-ZERO = TRUE
*  For AND logic,
*	All conditions must be true for a TRUE result
*	If any condition is flase the result is FALSE
*  For OR logic,
*	If any condition is true the result is TRUE
*	All conditions must be false for a FALSE result
*/

}

