/* Program 4.3 : Loan Repayment */

#include <stdio.h>

void main(void)
{
  float debt, MonthlyRate, payment;
  int month;

  printf("\nDebt:");  scanf("%f", &debt);
  printf("\nMonthly rate of interest:");  scanf("%f", &MonthlyRate);
  printf("\nMonthly payment:");  scanf("%f", &payment);

  for (month = 1; month <= 12; month++)
  {
    debt = debt + debt*MonthlyRate/100 - payment;
    printf("Debt after next payment is %.2f\n", debt);
  }
}

