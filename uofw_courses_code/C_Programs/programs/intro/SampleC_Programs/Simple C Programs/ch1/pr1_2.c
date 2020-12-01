/* Program 1.2 : Room Size 2 */

#include <stdio.h>

void main (void)
{
  int length, width, height;

  printf("Type in length, width and height:\n");
  scanf("%i%i%i", &length, &width, &height);
  printf("Your room needs %i square metres of carpet\n",
                                            length*width);
  printf("and %i square metres of wallpaper.\n",
                                 2*(length+width)*height);

}

