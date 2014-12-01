/* Program 3.8 : Loan */

#include <stdio.h>

void main(void)
{
  int loan;
  float InterestRate;

  printf("Amount of loan?");  scanf("%i", &loan);

  switch (loan / 1000)
  {
    case 0:  InterestRate = 10;    break;
    case 1:  InterestRate = 11;    break;
    case 2:  InterestRate = 11.5;  break;
    case 3:  InterestRate = 11.75; break;
    case 4:  InterestRate = 12;    break;
    default: InterestRate = 13;
  }
  printf("Interest to pay: %.2f\n", loan * InterestRate/100);

}

