/* Program 5.6 : Loan Statement */

#include <stdio.h>

void main(void)
{
  float debt, MonthlyRate, payment;
  int month;

  printf("\nDebt:");  scanf("%f", &debt);
  printf("\nMonthly rate of interest:");  scanf("%f", &MonthlyRate);
  printf("\nMonthly payment:");  scanf("%f", &payment);

  printf("\nMonth   Outstanding debt\n\n");

  month = 0;

  while ( (month < 12) &&
          (debt + debt*MonthlyRate/100 >= payment) )
  {
    month++;
    debt += debt*MonthlyRate/100 - payment;
    printf("%5i%19.2f\n", month, debt);
  }

  printf("\n");
  if (month < 12)
  {
    printf("Debt cleared within a year.\n");
    printf("Final payment: ");
    printf("%.2f\n", debt + debt*MonthlyRate/100);
  }
}

