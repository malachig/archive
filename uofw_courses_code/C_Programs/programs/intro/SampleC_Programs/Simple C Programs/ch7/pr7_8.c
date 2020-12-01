/* Program 7.8 : Add */

#include <stdio.h>

float max(float first, float second);

void main(void)
{ float a, b, p, q, x, y;

  scanf("%f%f%f%f%f%f", &a, &b, &p, &q, &x, &y);
  printf("%f\n", max(a, b) + max(p, q) + max(x, y));

}

float max(float first, float second)
{
  if (first > second)
    return first;
  else
    return second;
}

