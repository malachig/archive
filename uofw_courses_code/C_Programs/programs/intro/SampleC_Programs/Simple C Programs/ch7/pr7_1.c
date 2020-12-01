/* Program 7.1 : Exam Marks 2 */
#include <stdio.h>

void ProcessCandidate(void);

int candidate;

void main(void)
{
  for (candidate=1; candidate<=25; candidate++)
     ProcessCandidate();

}

void ProcessCandidate(void)
{ int exam, mark, total;
  float average;

    total =0;
    for (exam=1; exam<=5; exam++)
    {
      scanf("%i", &mark);
      total += mark;
    }
    average = total / 5.0;
    printf("Candidate %i", candidate);
    printf(" - Average mark: %.1f ", average);
    if (average >= 50)
      printf(" [Pass]\n");
    else
      printf(" [Fail]\n");
}

