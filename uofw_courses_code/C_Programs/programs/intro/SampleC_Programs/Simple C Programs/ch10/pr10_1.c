/* Program 10.1 : Introducing structs */
#include <stdio.h>
#include <string.h>


typedef struct {
  char title[51];
  int year;
  float price;
} books;

void main (void) {

  books zen, poker;

  strcpy(zen.title,"Zen and the Art of Motorcycle maintainance");
  zen.year = 1974;
  zen.price = 4.99;

  strcpy(poker.title,"The Education of a Poker Player");
  poker.year = 1957;
  poker.price = 5.95;

  printf("%s was published in %i and costs %.2f\n",
                            zen.title, zen.year, zen.price);
  printf("%s was published in %i and costs %.2f\n",
                            poker.title, poker.year, poker.price);
}

