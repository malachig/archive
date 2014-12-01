/* Program 10.4 : Passing a struct to a function */

#include <stdio.h>

typedef struct {
  int pages;
  float price;
} books;

void ShowPrice ( books SomeBook );

void main (void) 
{
  books zen;
  zen.pages = 400;
  zen.price = 4.99;
  printf("Zen costs ");
  ShowPrice(zen);
}

void ShowPrice ( books SomeBook ) 
{
  printf("%.2f\n", SomeBook.price);
}

