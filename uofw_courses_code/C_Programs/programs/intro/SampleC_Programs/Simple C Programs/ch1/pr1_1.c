/* Program 1.1 : Room Size */

#include <stdio.h>

void main (void)
{
  int length, width, height;

  scanf("%i%i%i", &length, &width, &height);
  printf("Your room needs %i square metres of carpet\n",
                                            length*width);
  printf("and %i square metres of wallpaper.\n",
                                 2*(length+width)*height);

}

