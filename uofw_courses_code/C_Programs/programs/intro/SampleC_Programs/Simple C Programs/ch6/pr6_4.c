/* Program 6.4 : Coins */

#include <stdio.h>

void main(void)
{
  int weight, total;
  total = 0;

  do {
    scanf("%i", &weight);
    switch (weight)
    {
      case 35: total += 50; break;
      case 16: total += 10; break;
      case 19: total += 20; break;
      case  9: total += 5;  break;
      case  7: total += 2;  break;
      case  3: total++;     break;
      default: printf("Coin rejected\n");
    }
  } while (total < 123);

  printf("Coins accepted.\n");
  if (total > 123)
    printf("Change due: %i\n", total - 123);

}

