/* Program 8.2 : Anagrams */

#include <stdio.h>

void main(void)
{
  char first, second, third;

  do {
    first = getchar();
  } while (first == ' ');

  second = getchar(); third  = getchar();

  printf("%c%c%c  ", first, second, third);
  printf("%c%c%c\n", first, third, second);

  printf("%c%c%c  ", second, first, third);
  printf("%c%c%c\n", second, third, first);

  printf("%c%c%c  ", third, first, second);
  printf("%c%c%c\n", third, second, first);

}

