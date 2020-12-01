/* Program 10.7 : Sort */

#include <stdio.h>

typedef struct {
  char title [51];
  int year;
  float price;
} books;

FILE *BookFile;
books book[21];
int LastBook;

void LoadBookFile    (void);
void LoadNextBook (books *onebook);
void DisplayAllBooks (void);
void DisplayOneBook  (books b);
void SwapBooks (int FirstIndex, int SecondIndex);

void main(void)
{
  int NumberToSort, LargestKey, LargestKeyWasAt, index;

  LoadBookFile();
  printf("\nBefore sort:\n");
  DisplayAllBooks();

  for (NumberToSort = LastBook; NumberToSort >= 2; NumberToSort--)
  {
    LargestKey = book[1].year;
    LargestKeyWasAt = 1;
    for (index = 2; index <= NumberToSort; index++)
    {
      if (book[index].year > LargestKey)
      {
        LargestKey = book[index].year;
        LargestKeyWasAt = index;
      }
    }
    if (LargestKeyWasAt != NumberToSort)
      SwapBooks (LargestKeyWasAt, NumberToSort);
  }

  printf("\nAfter sort:\n");
  DisplayAllBooks();
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

void DisplayAllBooks (void)
{
  int i;
  for (i=1; i<=LastBook; i++)
  {
    printf ("%2i ", i);
    DisplayOneBook (book[i]);
  }
}

void DisplayOneBook (books b)
{
  printf ("%s %4i %7.2f\n", b.title, b.year, b.price );
}

void SwapBooks (int FirstIndex, int SecondIndex)
{
  books temp;
  temp = book[FirstIndex];
  book[FirstIndex] = book[SecondIndex];
  book[SecondIndex] = temp;
}

