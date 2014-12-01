/* Program 5.10 : Election 2 */

#include <stdio.h>

void main(void)
{
  int Party1Next, Party2Next,
      Party1Overall, Party2Overall;

  Party1Overall = 0;
  scanf("%i", &Party1Next);
  while (Party1Next >= 0)
  {
    Party1Overall += Party1Next;
    scanf("%i", &Party1Next);
  }

  Party2Overall = 0;
  scanf("%i", &Party2Next);
  while (Party2Next >= 0)
  {
    Party2Overall += Party2Next;
    scanf("%i", &Party2Next);
  }

  printf("Party 1: %i\n", Party1Overall);
  printf("Party 2: %i\n", Party2Overall);
    
}

