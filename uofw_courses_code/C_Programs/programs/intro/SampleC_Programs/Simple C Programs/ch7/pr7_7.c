/* Program 7.7 : Change  - Version 2*/

#include <stdio.h>

void HowMany(int CoinValue, int *AmountLeft);

void main (void)
{
  int price, change;

  printf("\nPrice in pence(<100):");
  scanf("%i", &price);
  change = 100 - price;

  printf("Change due is:\n");
  HowMany(50, &change);
  HowMany(20, &change);
  HowMany(10, &change);
  HowMany(5, &change);
  HowMany(2, &change);
  HowMany(1, &change);
  printf("\n");

}

void HowMany(int CoinValue, int *AmountLeft)
{  int coins;
   coins = *AmountLeft/CoinValue;
   *AmountLeft = *AmountLeft % CoinValue;
   printf("%i %is, ", coins, CoinValue);
}

