/* Program 3.6 : Jury Service */

#include <stdio.h>

void main(void)
{
  int age;
  printf("What is your age?");
  scanf("%i", &age);

  if ( (age>=18) && (age<=65) )
    printf("Eligible for jury service.\n");
  else
    printf("Not eligible for jury service.\n");

}

