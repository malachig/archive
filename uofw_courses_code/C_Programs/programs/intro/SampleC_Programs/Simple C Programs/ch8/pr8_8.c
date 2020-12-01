/* Program 8.8 : Count Four-Letter Words */

#include <stdio.h>
#include <ctype.h>


void FindStartOfWord(void);
int LettersInWord(void);

char nextch;

void main(void)
{
  int letters, words, fourwords;
  words = 0; fourwords = 0;
  do
  {
    FindStartOfWord();
    letters = LettersInWord();
    words++;
    if (letters == 4) fourwords++;
  } while (nextch != '.');

  printf("Words in sentence: %i\n", words);
  printf("4-letter words: %i\n", fourwords);
}


void FindStartOfWord(void)
{
  do
  {
    nextch = getchar();
  } while (!isalpha(nextch));
}


int LettersInWord(void)
{
  int count;
  count = 0;

  do
  {
    count++;
    nextch = getchar();
  } while (isalpha(nextch));

  return count;
}

