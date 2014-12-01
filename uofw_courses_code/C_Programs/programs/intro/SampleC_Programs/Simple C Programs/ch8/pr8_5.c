/* Program 8.5 : Classify Character */

#include <stdio.h>
#include <ctype.h>

void main(void)
{
  char character;
  printf("Type character for analysis:");
  character = getchar();

  if (isalpha(character))
    printf("That character is alphabetic.\n");
  else
    if (isdigit(character))
      printf("That character is numeric.\n");
    else
      printf("That character is not alphanumeric.\n");

}

