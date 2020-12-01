/* Program 8.11 : Add Up */

#include <stdio.h>
#include <ctype.h>

void main(void)
{
  int HowMany, count, next, total;
  char c;

  /* 'First' section */
  do
    c = getchar();
  while (!isdigit(c));
  ungetc(c, stdin);
  scanf("%i", &HowMany);

  /* 'Second' section */
  total = 0;
  for (count = 1; count <= HowMany; count++)
  {
    do
      c = getchar();
    while (!isdigit(c));
    ungetc(c, stdin);
    scanf("%i", &next);
    total += next;
  }

  /* 'Third' section */
  do c = getchar(); while (c != '?');

  printf("Total is: %i\n", total);

}

