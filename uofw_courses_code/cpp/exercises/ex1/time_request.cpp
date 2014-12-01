#include <iostream>
#include <string>
using namespace std;

int main()
{
	int day, month, year;
	string weekday;  //string data type allowed by string stream above.

	cout << "Enter the day > ";
	cin >> day;
	
	cout << "Enter the month > ";
	cin >> month;

	cout << "Enter the year > ";
	cin >> year;

	cout << "Enter the day of the week > ";
	cin >> weekday;

	cout << endl << "Printing in Reverse Order:" << endl;
	cout << weekday << endl << year << endl << month << endl << day << endl;

	return 0;
}
