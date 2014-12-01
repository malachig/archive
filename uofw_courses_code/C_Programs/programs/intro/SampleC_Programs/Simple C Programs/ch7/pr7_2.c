/* Program 7.2 : Exam Marks 3 */
#include <stdio.h>

void ProcessCandidate(int cand);

void main(void)
{
  int candidate;
  for (candidate=1; candidate<=25; candidate++)
     ProcessCandidate(candidate);

}

void ProcessCandidate(int cand)
{ int exam, mark, total;
  float average;

    total =0;
    for (exam=1; exam<=5; exam++)
    {
      scanf("%i", &mark);
      total += mark;
    }
    average = total / 5.0;
    printf("Candidate %i", cand);
    printf(" - Average mark: %.1f ", average);
    if (average >= 50)
      printf(" [Pass]\n");
    else
      printf(" [Fail]\n");
}

