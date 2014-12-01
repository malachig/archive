/* Program 10.2 : Latest books */

#include <stdio.h>

typedef struct {
  char title [51];
  int year;
  float price;
} books;

FILE *BookFile;
books book[21];
int LastBook;

void LoadBookFile (void);

void main (void)
{
  int index, LatestYear;

  LoadBookFile();
  LatestYear = book[1].year;
  for (index = 2; index <= LastBook; index++)
    if (book[index].year > LatestYear)
      LatestYear = book[index].year;

  for (index = 1; index <= LastBook; index++)
  {
    if (book[index].year == LatestYear)
      printf ("*");
    else
      printf (" ");

    printf ("%3i. ", index);

    printf ("%s %4i %7.2f\n", book[index].title,
       book[index].year, book[index].price );
  }

}

void LoadBookFile (void)
{
  char NextChar;
  int NextFree, i, year;
  float price;

  BookFile = fopen("BOOKS.DAT", "r");
  NextFree = 1;
  NextChar = getc(BookFile);

  while (!feof(BookFile) && NextFree <= 20)
  {
    book[NextFree].title[0] = NextChar;
    for (i=1; i<=49; i++)
      book[NextFree].title [i] = fgetc(BookFile);
    book[NextFree].title [50] = '\0';
    fscanf ( BookFile, "%i %f", &year, &price); 
    book[NextFree].year  = year;
    book[NextFree].price = price;
    NextFree++;
    do {
      NextChar = getc(BookFile);
    } while ( NextChar != '\n' && !feof(BookFile) );
    NextChar = getc(BookFile);
  }
  LastBook = NextFree - 1;
  fclose(BookFile);
}

