/* Program 3.7 : Insurance Policy */

#include <stdio.h>

void main(void)
{
  int age, cc, convictions;
  int over21, LargeCar, RiskDriver;

  printf("Age?"); scanf("%i", &age);
  printf("cc?"); scanf("%i", &cc);
  printf("How many convictions?"); scanf("%i", &convictions);

  over21 = (age>=21);
  LargeCar = (cc>=2000);
  RiskDriver = (convictions>=3);

  if (over21 && LargeCar && RiskDriver)
    printf("Policy loaded by 45 percent.\n");

  if (over21 && LargeCar && !RiskDriver)
    printf("Policy loaded by 15 percent.\n");

  if (over21 && !LargeCar && RiskDriver)
    printf("Policy loaded by 30 percent.\n");

  if (over21 && !LargeCar && !RiskDriver)
    printf("No loading.\n");

  if (!over21 && LargeCar && RiskDriver)
    printf("No policy to be issued.\n");

  if (!over21 && LargeCar && !RiskDriver)
    printf("Policy loaded by 60 percent.\n");

  if (!over21 && !LargeCar && RiskDriver)
    printf("Policy loaded by 50 percent.\n");

  if (!over21 && !LargeCar && !RiskDriver)
    printf("Policy loaded by 10 percent.\n");

}

