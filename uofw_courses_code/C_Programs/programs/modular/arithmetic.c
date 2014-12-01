/* Author: Malachi Griffith
*  Date: Nov. 7 2002
*  Purpose: Adds, subtracts, multiplies and divides common fractions, 
*  displaying results in reduced form.
*/

#include <stdio.h>
#include <stdlib.h>  /* provides function abs */
#include <math.h>

/* Function prototypes */
void scan_fraction (int *nump, int *denomp);
char get_operator(void);
void add_fractions(int n1, int d1, int n2, int d2, int *n_ansp, int *d_ansp);
void multiply_fractions(int n1, int d1, int n2, int d2,
			int *n_ansp, int *d_ansp);
int find_gcd (int n1, int n2);
void reduce_fraction(int *nump, int *denomp);
void print_fraction(int num, int denom);

int
main(void)
{
	int n1, d1;	/* numerator, denominator of first fraction */
	int n2, d2;	/* numerator, denominator of second fraction */
	char op;	/* arithmetic operator + - * or /  */
	char again;	/* y or n depending on user's desire to continue */
	int n_ans,	/* numerator, denominator of answer */
	    d_ans;

	/* While the user wants to continue, gets and solves arithmetic
	   problems with common fractions. */

	do
	{
       	   /* Gets a fraction problem */
	   scan_fraction(&n1, &d1);
	   op = get_operator();
	   scan_fraction(&n2, &d2);

	   /* Computes the result */
	   switch (op)
	   {
	  	case '+':
			add_fractions(n1, d1, n2, d2, &n_ans, &d_ans);
			break;

		case '-':
			add_fractions(n1, d1, -n2, d2, &n_ans, &d_ans);
			break;

		case '*':
			multiply_fractions(n1, d1, n2, d2, &n_ans, &d_ans);
			break;

		case '/':
			multiply_fractions(n1, d1, d2, n2, &n_ans, &d_ans);
	   }
	   reduce_fraction(&n_ans, &d_ans);

	   /* Displays problem and result */
	   printf("\n");
	   print_fraction(n1, d1);
	   printf(" %c ", op);
	   print_fraction(n2, d2);
	   printf(" = ");
	   print_fraction(n_ans, d_ans);

	   /* Asks user about doing another problem */
	   printf("\nDo another problem? (y/n) > ");
	   scanf(" %c", &again);
	}
	while (again == 'Y' || again == 'y');
	return(0);
}

/*  Function scan_fraction */
/* 
*   Gets and returns a valid fraction as its result.
*   A valid fraction is of this form:  integer/positive integer.
*/
void
scan_fraction(int *nump, int *denomp)
{
	char slash;	/* character between numerator and denominator */
	int status;	/* status code returned by scanf indicating number
			 * of valid values obtained */
	int error;	/* flag indicating presence of an error */
	char discard;	/* unprocessed character from input file */

	do
	{
	   /*No errors detected yet */
	   error = 0;
	   
	   /* Get a fraction from the user */
	   printf("Enter a common fraction as two integers separated ");
	   printf("by a slash > ");
	   status = scanf("%d %c%d", nump, &slash, denomp);

	   /* Validate the fraction */
	   if (status < 3)
	   {
		error = 1;
		printf("Invalid-please read directions carefully\n");
	   }
	   else if (slash != '/')
		{
		error = 1;
		printf("Invalid-seperate numerator and denominator");
		printf(" by a slash (/)\n");
		}
	   else if (*denomp <= 0)
	 	{
		error = 1;
		printf("Invalid-denominator must be positive\n");
		}
	   /* Discard extra input characters */
	   do 
	   {
		scanf("%c", &discard);
	   }
	   while (discard != '\n');
	}
	while (error);
}

/* 
*  Gets and returns a valid arithmetic operator.  Skips over newline 
*  characters and permits re-entry of operator in case of error.
*/
char
get_operator(void)
{
	char op;

	printf("Enter an arithmetic operator (+, -, *, /)\n >");
	for (scanf("%c", &op);
	     op != '+' && op != '-' && op != '*' && op != '/';
	     scanf("%c", &op))
	{
		if (op != '\n')
			printf("%c invalid, reenter operator(+, -, *, /)\n >",
			       op);
	}
	return(op);
}
 
/*
*  Adds fractions represented by pair of integers.
*  Pre: n1, d1, n2, d2 are defined;
*       n_ansp and d_ansp are addresses of type int variables.
*  Post: sum of n1/d1 and n2/d2 is sored in variables pointed
*        to by n_ansp and d_ansp.  Result is not reduced.
*/
void
add_fractions (int n1, int d1, int n2, int d2, int *n_ansp, int *d_ansp)
{
	int denom, 	/*common denominator used for sum (may not be least)*/
	numer,		/* numerator of sum */
	sign_factor;	/* -1 for a negative, 1 otherwise. */

	/* Finds a common denominator */
	denom = d1 * d2;

	/* Computes numerator */
	numer = n1 * d2 + n2 * d1;

	/* Adjusts sign (at most, numerator should be negative) */
	if (numer * denom >= 0)
		sign_factor = 1;
	else
 		sign_factor = -1;

	numer = sign_factor * abs(numer);
	denom =abs(denom);

	/* Returns result */
	
	*n_ansp = numer;
	*d_ansp = denom;
}

/*
* Multiplies fractions represented by pairs of integers.
* Pre: n1, d1, n2, d2 are defined;
*      n_ansp and d_ansp are addresses of type int variables.
* Post: product of n1/d1 and n2/d2 is stored in variables pointed to 
*       by n_ansp and d_ansp.  Result is not reduced.
*/
void
multiply_fractions(int n1, int d1, int n2, int d2, int *n_ansp, int *d_ansp)
{
	int numer, denom;

	/* Displays trace message */
	printf("\nEntering multiply fractions with\n");
	printf("n1 = %d, d1 = %d, n2 = %d, d2 = %d\n", n1, d1, n2, d2);

	numer = n2 * n2;
	denom = d1 * d2;	

	/* Defines output arguments */
	
	*n_ansp = numer;
	*d_ansp = denom;
}

/*
*  Finds the greatest common divisor of two integers 
*/
int 
find_gcd (int n1, int n2)
{
	int gcd;
	int p, q, r;

	q = abs(n1);
	p = abs(n2);
	
	r = q % p;
	
	while (r != 0)	
	{
		q = p;
		p = r;
		r = q % p;
	}
	
	gcd = p;

	/* Displays trace messages */
	printf("\nEntering find_gcd with n1 = %d, n2 = %d\n", n1, n2);

	/* Displays exit trace message */
	printf("find_gcd returning %d\n", gcd);
	return (gcd);
}

/*
*  Reduces a fraction by diving its numerator and denominator by their
*  greatest common divisor
*/
void
reduce_fraction (int *nump, int *denomp)
{
	int gcd;  /* greatest common divisor of numerator & denominator */

	gcd = find_gcd (*nump, *denomp);

	*nump = *nump / gcd;
	*denomp = *denomp / gcd;
}

/*
*  Displays pair of integers as a fraction.
*/
void
print_fraction (int num, int denom)  /* input - numerator */
{
	printf("%d/%d", num, denom);
}


