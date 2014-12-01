/* Program 5.5 : Arithmetic Test */

#include <stdio.h>
#define A 16
#define B 25

void main(void)
{
  int attempts, answer;
  printf("Let's test your arithmetic.\n");
  printf("Type the answer to the following sum.\n");
  printf("%i + %i = ", A, B);
  scanf("%i", &answer); attempts = 1;
  while (answer != A+B)
  {
    printf("\nWrong - try again.\n");
    printf("%i + %i = ", A, B);
    scanf("%i", &answer); attempts++;
  }
  printf("\n");
  if (attempts == 1)
    printf("Very good - you got it in one!\n");
  else printf("You got it at last!\n");

  printf("Bye for now!\n");
}

