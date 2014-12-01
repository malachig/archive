/* Program 8.9 : Code */

#include <stdio.h>
#include <ctype.h>

char coded(char ch);

void main(void)
{
  char character;
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
{
    if (isalpha(ch))
      if (ch == 'z')
        return 'a';
      else
        return ch+1;
    else
      return ch;
}

