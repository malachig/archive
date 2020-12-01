/* Program 8.10 : Code Shift */

#include <stdio.h>
#include <ctype.h>
#define SHIFT 5

char coded(char ch);

void main(void)
{
  char character; int ordinal;
  FILE *text;
  text = fopen("TEXT.DAT","r");

  character = getc(text);
  while (!feof(text))
  {
    putchar(coded(character));
    character = getc(text);
  } 

  fclose(text);
  }

char coded(char ch)
{ char newch;
    if (isalpha(ch))
    {
      newch = ch + SHIFT;
      if (newch > 'z') newch -= 26;
    }
    else newch = ch;
    return newch;
}

