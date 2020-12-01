/* Program 7.4 : Compare */

#include <stdio.h>

void CommaWrite(long int n);
void ZeroWrite(int m);

void main(void)
{ long int first, second;

  scanf("%li%li", &first, &second);
  if (first > second) 
  { CommaWrite(first); CommaWrite(second); }
  else
  { CommaWrite(second); CommaWrite(first); }
}

void CommaWrite(long int n)
{
  int millions, thousands, units;

  millions = (n / 1000000);
  thousands = (n % 1000000) / 1000;
  units = (n % 1000);

  if (millions > 0)
  {
    printf("%i,", millions);
    ZeroWrite(thousands); putchar(',');
    ZeroWrite(units);
  }
  else
  {
    if (thousands > 0)
    {
      printf("%i,", thousands);
      ZeroWrite(units);
    }
    else printf("%i", units);
  }
  printf("\n");
}

void ZeroWrite(int m)
{
  if (m < 100) printf("0");
  if (m < 10)  printf("0");
  printf("%i", m);
}

