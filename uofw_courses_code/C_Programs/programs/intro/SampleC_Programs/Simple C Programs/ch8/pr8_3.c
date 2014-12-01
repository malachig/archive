/* Program 8.3 : Framed Triangle 2 */

#include <stdio.h>
void PrintRow(int row, int rows);
void PrintChars(char ch,  int n);

void main(void)
{
  int rows, WidthOfPicture, NextRow;
  printf("\nheight of triangle:");  scanf("%i", &rows);
  WidthOfPicture = 2*rows + 1;

  PrintChars('.',WidthOfPicture);
  printf("\n");

  for(NextRow=1; NextRow<=rows; NextRow++)
    PrintRow(NextRow, rows);

  PrintChars('.',WidthOfPicture);
  printf("\n");

}


void PrintRow(int row, int rows)
{
  PrintChars('.', rows + 1 - row);
  PrintChars('*', 2*row - 1);
  PrintChars('.', rows + 1 - row);
  printf("\n");
}


void PrintChars(char ch,  int n)
{
  int count;
  for (count=1; count<=n; count++) putchar(ch);
}

