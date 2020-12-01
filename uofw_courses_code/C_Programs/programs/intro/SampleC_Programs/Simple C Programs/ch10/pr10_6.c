/* Program 10.6 : Book search */

#include <stdio.h>

#define TRUE 1
#define FALSE 0

typedef struct {
  char title [51];
  int year;
  float price;
} books;

FILE *BookFile;
books book[21];
int LastBook;

void LoadBookFile (void);
void LoadNextBook (books *onebook);
void DisplayMenu (void);
void SearchByYear (void);
void SearchByPrice (void);
void DisplayOneBook (books b);

void main (void)
{
  int command;
  LoadBookFile();
  do {
    DisplayMenu();
    scanf ("%i", &command);
    if (command == 1) SearchByYear();
    if (command == 2) SearchByPrice();
  } while (command != 3);
}

void LoadBookFile (void)
{
  int NextFree;
  books BookBuffer;
  BookFile = fopen("BOOKS.DAT", "r");
  NextFree = 1;
  LoadNextBook (&BookBuffer);
  while (!feof(BookFile) && NextFree <= 20)
  {
    book[NextFree] = BookBuffer;
    NextFree++;
    LoadNextBook (&BookBuffer);
  }
  LastBook = NextFree - 1;
  fclose(BookFile);
}

void LoadNextBook (books *onebook)
{
  int i, ch;
  int year;
  float price;

  for (i=0; i<=49; i++)
    (*onebook).title [i] = getc(BookFile);
  (*onebook).title [50] = '\0';

  fscanf (BookFile, "%i %f", &year, &price);

  (*onebook).year = year;
  (*onebook).price = price;

  do {
    ch = getc(BookFile);
  } while ( ch != '\n' && !feof(BookFile) );
}

void DisplayMenu (void)
{
  printf ("\n");
  printf ("1: Search by year\n");
  printf ("2: Search by price\n");
  printf ("3: Quit\n");
  printf ("Enter selection: ");
}

void SearchByYear (void)
{
  int SearchYear, SearchIndex, NoBooks;
  printf ("Year of publication: "); scanf  ("%i", &SearchYear);
  NoBooks = TRUE;

  for (SearchIndex=1; SearchIndex<=LastBook; SearchIndex++)
  {
    if (book[SearchIndex].year == SearchYear)
    {
      DisplayOneBook (book[SearchIndex]);
      NoBooks = FALSE;
    }
  }

  if (NoBooks)
    printf ("There are no books published that year.\n");
}

void SearchByPrice (void)
{
  int SearchIndex, NoBooks;
  float UpperLimit;

  printf ("Upper limit for price: "); scanf  ("%f", &UpperLimit);
  NoBooks = TRUE;

  for (SearchIndex=1; SearchIndex<=LastBook; SearchIndex++)
  {
    if (book[SearchIndex].price <= UpperLimit)
    {
      DisplayOneBook (book[SearchIndex]);
      NoBooks = FALSE;
    }
  }

  if (NoBooks)
    printf ("There are no books in your price-range.\n");
}

void DisplayOneBook (books b)
{
  printf ("%s %4i %3.2f\n", b.title, b.year, b.price );
}

