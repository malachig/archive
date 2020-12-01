/* Program 3.3 : Loan Repayment */

#include <stdio.h>
#include <math.h>

void main (void)
{
  float debt, MonthlyRate, payment;

  printf("debt?");  scanf("%f", &debt);
  printf("monthly rate?");  scanf("%f", &MonthlyRate);
  printf("repayment?");  scanf("%f", &payment);
  debt = debt + debt*MonthlyRate/100 - payment;

  if (debt < 0)
    printf("Debt cleared. Refund = %.2f\n", -debt);
  else
    printf("Debt after this payment is %.2f\n", debt);

}

