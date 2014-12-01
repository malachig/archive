/* Program 1.5 : ElectricityBill */

#include <stdio.h>

void main (void)
{
  int present, previous;

  printf("Type present and previous meter readings:\n");
  scanf("%i%i", &present, &previous);

  printf("*****************************\n");
  printf("Present meter reading    %i\n", present);
  printf("Previous meter reading   %i\n", previous);
  printf("Units used                %i\n", present - previous);
  printf("Rate per unit            6.73 pence\n");
  printf("Standing charge         £7.78\n\n");
  printf("The sum due is         £");
  printf("%f\n", (present-previous) * 6.73/100 + 7.78);
  printf("*****************************\n");

}

