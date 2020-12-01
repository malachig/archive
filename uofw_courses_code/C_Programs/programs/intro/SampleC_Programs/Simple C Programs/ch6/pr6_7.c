/* Program 6.7 : Sales Analysis */
#include <stdio.h>

void main(void)
{
  int dept, LargeSales;
  float NextItemPrice;

  for (dept = 1; dept <= 63; dept++)
  { LargeSales = 0;
    scanf("%f", &NextItemPrice);

    while (NextItemPrice >= 0.0)
    {
      if (NextItemPrice >= 10.0)
        LargeSales++;
      scanf("%f", &NextItemPrice);
    }

    printf("Department %i has made ", dept);
    printf("%i large sales.\n", LargeSales);
  }

}

