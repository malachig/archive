/* Program 6.6 : Framed Triangle */

#include <stdio.h>

void main(void)
{
  int rows, WidthOfPicture, count,
      row, stars, dots;

  printf("\nheight of triangle:");  scanf("%i", &rows);

  WidthOfPicture = 2*rows + 1;

  printf("\n");

  for (count=1; count<=WidthOfPicture; count++)
    printf(".");
  printf("\n");

  for(row=1; row<=rows; row++)
  { stars = 2*row - 1;
    dots  = rows + 1 - row;
    for (count=1; count<=dots; count++) printf(".");
    for (count=1; count<=stars; count++) printf("*");
    for (count=1; count<=dots; count++) printf(".");
    printf("\n");
  }

  for (count=1; count<=WidthOfPicture; count++)
    printf(".");
  printf("\n");

}

