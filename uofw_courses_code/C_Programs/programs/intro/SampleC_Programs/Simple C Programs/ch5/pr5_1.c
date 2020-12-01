/* Program 5.1 : Complete Statement */

#include <stdio.h>

void main(void)
{
  float debt, MonthlyRate, payment;
  printf("\nDebt:");  scanf("%f", &debt);
  printf("\nMonthly rate of interest:");  scanf("%f", &MonthlyRate);
  printf("\nMonthly payment:");  scanf("%f", &payment);

  while (debt + debt*MonthlyRate/100 >= payment)
  {
    debt += debt*MonthlyRate/100 - payment;
    printf("Debt after next payment is %.2f\n", debt);
  }

  printf("\nFinal payment required will be %.2f\n",
                         debt + debt*MonthlyRate/100);

}

