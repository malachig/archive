/* Program 1.3 : LoanRepayment */

#include <stdio.h>

void main (void)
{
  float debt, MonthlyRate, payment;

  printf("debt?");  scanf("%f", &debt);
  printf("monthly rate?");  scanf("%f", &MonthlyRate);
  printf("repayment?");  scanf("%f", &payment);
  printf("Debt after next payment is ");
  printf("%f\n", debt+debt*MonthlyRate/100 - payment);

}

