/* Program 1.7 : ElectricityBill2 */

#include <stdio.h>
#define UNITRATE 6.73
#define STANDCHARGE 7.78
#define ASTERISKS "*****************************\n"

void main (void)
{
  int present, previous;

  printf("Type present and previous meter readings:\n");
  scanf("%i%i", &present, &previous);
  printf(ASTERISKS);
  printf("Present meter reading %7i\n", present);
  printf("Previous meter reading%7i\n", previous);
  printf("Units used            %7i\n", present - previous);
  printf("Rate per unit         %7.2f pence\n", UNITRATE);
  printf("Standing charge      £%7.2f\n\n", STANDCHARGE);
  printf("The sum due is       £%7.2f\n",
          (present-previous) * UNITRATE/100 + STANDCHARGE );
  printf(ASTERISKS);

}

