/* Program 8.4 : Character count */

#include <stdio.h>
#define GIVENCH 'e'

void main(void)
{
  int occurrences; char nextch;
  FILE *text;
  text = fopen("TEXT.DAT","r");

  occurrences = 0;

  nextch = getc(text);
  while (!feof(text))
  {
    if (nextch == GIVENCH) occurrences++;
    nextch = getc(text);
  } 

  fclose(text);

  printf("The character %c", GIVENCH);
  printf(" appears %i times.\n", occurrences);

}

