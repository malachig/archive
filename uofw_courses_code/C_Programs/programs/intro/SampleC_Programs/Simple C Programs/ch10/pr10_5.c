/* Program 10.5 : Price raise */

#include <stdio.h>

typedef struct {
  int pages;
  float price;
} books;

void RaisePrice ( books *SomeBook );

void main (void) {
  books zen;
  zen.pages = 400;
  zen.price = 4.99;
  printf ("Zen costs: %.2f\n", zen.price);
  RaisePrice(&zen);
  RaisePrice(&zen);
  printf ("Zen now costs: %.2f\n", zen.price);
}

void RaisePrice ( books *SomeBook )
{
  (*SomeBook).price = 1.1 * (*SomeBook).price;
}

