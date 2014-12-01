/* Program 10.3 : Swap Shoes */

typedef struct {
  char sex;
  int size;
} shoes;

void main (void)
{
  shoes MyPair, YourPair, temp;

  MyPair.sex = 'M';
  MyPair.size = 11;

  YourPair.sex = 'F';
  YourPair.size = 8;

  printf("My shoes are: %c %i.\n", MyPair.sex, MyPair.size);
  printf("Your shoes are: %c %i.\n", YourPair.sex, YourPair.size);

  temp = YourPair;
  YourPair = MyPair;
  MyPair = temp;

  printf("Now my shoes are: %c %i.\n", MyPair.sex, MyPair.size);
  printf("Now your shoes are: %c %i.\n", YourPair.sex, YourPair.size);

}

