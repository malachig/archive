/* Program 6.10 : Holiday */

#include <stdio.h>
#define LOWRATE  100
#define MIDRATE  200
#define PEAKRATE 300

void main(void)
{
  int day1, month1, day2, month2;
  int duration, TillEndOfMonth, cost;

  scanf("%i%i", &day1, &month1);
  scanf("%i%i", &day2, &month2);

  if ( (day1>=1) && (day1<=31) &&
       (day2>=1) && (day2<=31) &&
       (month1>=1) && (month1<=12) &&
       (month2>=1) && (month2<=12) &&
       (month1 <= month2) )

  {
    if (month1 == month2)
      duration = day2 - day1 + 1;
    else
    {
      switch(month1)
      {
        case 2  :
          TillEndOfMonth = 28 - day1; break;
        case 9  : case 4: case 6: case 11:
          TillEndOfMonth = 30 - day1; break;
        default :
          TillEndOfMonth = 31 - day1;
      }
      duration = day2 + TillEndOfMonth + 1;
    }

    switch(month2) {
      case 1: case 2: case 11: case 12:
        cost = duration * LOWRATE; break;
      case 3: case 4: case 5: case 10:
        cost = duration * MIDRATE; break;
      case 6: case 7: case 8: case 9:
        cost = duration * PEAKRATE; break;
    }
    printf("Cost of holiday is $%i\n", cost);
  }
  else printf("You have typed erroneous dates.\n");

}

