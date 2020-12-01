/* Program 4.4 : Loan Statement */

#include <stdio.h>

void main(void)
{
  float debt, MonthlyRate, payment;
  int month;

  printf("\nDebt:");  scanf("%f", &debt);
  printf("\nMonthly rate of interest:");  scanf("%f", &MonthlyRate);
  printf("\nMonthly payment:");  scanf("%f", &payment);

  printf("\nmonth   outstanding debt\n");

  for (month = 1; month <= 12; month++)
  {
    debt = debt + debt*MonthlyRate/100 - payment;
    printf("%5i%19.2f\n", month, debt);
  }

}

