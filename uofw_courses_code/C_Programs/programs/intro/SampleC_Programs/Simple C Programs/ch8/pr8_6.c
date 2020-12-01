/* Program 8.6 : Word Length */

#include <stdio.h>
#include <ctype.h>

void main(void)
{
  char NextChar; int LetterCount;

  do {
    NextChar = getchar();
  } while (!isalpha(NextChar));

  /* At this point, NextChar will contain */
  /* the first letter of the word.        */

  LetterCount = 0;
  do {
    LetterCount++;
    NextChar = getchar();
  } while (isalpha(NextChar));

  /* At this point, NextChar will contain */
  /* the first character after the word.  */

  printf("No. of letters: %i\n", LetterCount);

}

