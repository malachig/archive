/* Program 8.7 : Letter Frequency */
#include <stdio.h>
#include <ctype.h>
#define GIVEN_LETTER 'e'

void main(void)
{
  int LetterCount, GivenLetterCount;
  char character;
  FILE *text;
  text = fopen("TEXT.DAT","r");

  LetterCount = 0; GivenLetterCount = 0;
  character = getc(text);
  while (!feof(text))
  {
    if (isalpha(character))
    {
      LetterCount++;
      if (character == GIVEN_LETTER)
        GivenLetterCount++;
    }
    character = getc(text);
  } 

  fclose(text);

  printf("%f percent of the letters were %c.\n", 
              100.0 * GivenLetterCount/LetterCount, GIVEN_LETTER);

}

