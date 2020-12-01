/* Program 6.5 : Exam Marks */

#include <stdio.h>

void main(void)
{
  int candidate, exam, total, mark;
  float average;

  for (candidate=1; candidate<=25; candidate++)
  {
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

}

