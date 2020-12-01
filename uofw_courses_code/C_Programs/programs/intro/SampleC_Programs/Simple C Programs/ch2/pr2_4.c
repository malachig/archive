/* Program 2.4 : Change */

#include <stdio.h>

void main (void)
{
  int price, change,
      fifties, twenties, tens, fives, twos, ones;

  printf("\nPrice in pence(<100):");
  scanf("%i", &price);
  change = 100 - price;

  fifties = change / 50; change = change % 50;
  twenties = change / 20; change = change % 20;
  tens = change / 10; change = change % 10;
  fives  = change / 5;  change = change % 5;
  twos  = change / 2;  change = change % 2;
  ones  = change;

  printf("Change due is:\n");
  printf("%i 50s, %i 20s, %i 10s, %i 5s, %i 2s and %i 1s.\n",
        fifties, twenties, tens, fives, twos, ones );
}

