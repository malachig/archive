#include <iostream>
using namespace std;

#include "time.h"

main()
{
	Time t;		// Instantiate a Time object.
	Time t2;

	cout << "The initial military time is ";
	t.print_military();

	cout << "The initial standard time is ";
	t.print_standard();

	t.set_time(13, 27, 6);
	cout << endl << endl << "After setting the time to 13:27:06 " << endl;
	cout << "The military time is ";
	t.print_military();
	cout << "The standard time is ";
	t.print_standard();

	t2.set_time(99, 99, 99);
	cout << endl << endl << "After setting the time to 99:99:99 " << endl;
	cout << "The military time is ";
	t2.print_military();
	cout << "The standard time is ";
	t2.print_standard();

	return 0;
}
