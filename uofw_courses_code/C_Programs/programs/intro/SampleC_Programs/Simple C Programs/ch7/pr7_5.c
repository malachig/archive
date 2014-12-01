/* Program 7.5 : Election 3 */

#include <stdio.h>

void AddUpVotesFor(int *total);

void main(void)
{ int Party1Overall, Party2Overall;

  AddUpVotesFor(&Party1Overall);
  AddUpVotesFor(&Party2Overall);

  if (Party1Overall == Party2Overall)
    printf("A draw.\n");
  else if (Party1Overall > Party2Overall)
    printf("A win for party 1.\n");
  else
    printf("A win for party 2.\n");
}

void AddUpVotesFor(int *total)
{
  int next;
  *total = 0;
  scanf("%i", &next);
  do {
    *total += next;
    scanf("%i", &next);
  } while (next >= 0);
}

