/* Program 6.8 : Policy 2 */

#include <stdio.h>

#define P45 "Policy loaded by 45 percent.\n"
#define P15 "Policy loaded by 15 percent.\n"
#define P30 "Policy loaded by 30 percent.\n"
#define OK  "No loading.\n"
#define NO  "No policy to be issued.\n"
#define P60 "Policy loaded by 60 percent.\n"
#define P50 "Policy loaded by 50 percent.\n"
#define P10 "Policy loaded by 10 percent.\n"

void main(void)
{
  int over21, LargeCar, RiskDriver;
  int age, cc, convictions;

  scanf("%i%i%i", &age, &cc, &convictions);
  over21 = (age >= 21);
  LargeCar = (cc >= 2000);
  RiskDriver = (convictions >= 3);

  if (over21)
  {
    if (LargeCar) {
      if (RiskDriver) printf(P45);
      else printf(P15); }
    else {
      if (RiskDriver) printf(P30);
      else printf(OK); }
  }
  else
  {
    if (LargeCar) {
      if (RiskDriver) printf(NO);
      else printf(P60); }
    else {
      if (RiskDriver) printf(P50);
      else printf(P10); }
  }

}

