#include <iostream>
using namespace std;

#include "time.h"

// Default constructor.
Time::Time()
{
	// Set hour/min/sec to 0:0:0
	dm_hour = dm_minute = dm_second = 0;
}

// Method for setting time.
void Time::set_time(int hour, int minute, int second)
{
	// If time is not valid set to 0:0:0
	dm_hour = dm_minute = dm_second = 0;

	if (hour >= 0 || hour <= 23)
		dm_hour = hour;
	if (minute >= 0 || minute <= 59)
		dm_minute = minute;
	if (second >= 0 || second <= 59)
		dm_second = second;
}

// Copy constructor.
void Time::set_time(const Time &t)
{
	dm_hour = t.dm_hour;
	dm_minute = t.dm_minute;
	dm_second = t.dm_second;
}

// Display the time in military fashion.
void Time::print_military()
{
	if (dm_hour < 10)
		cout << "0";
	cout << dm_hour << ":";

	if (dm_minute < 10)
		cout << "0";
	cout << dm_minute << ":";

	if (dm_second < 10)
		cout << "0";
	cout << dm_second << endl;
}

// Display the time in standard fasion.
void Time::print_standard()
{
	if (dm_hour==0 || dm_hour==12)
		cout << dm_hour << ":";
	else 
		cout << dm_hour%12 << ":";

	if (dm_minute < 10)
		cout << "0";
	cout << dm_minute << ":";

	if (dm_second < 10)
		cout << "0";
	cout << dm_second;

	if (dm_hour < 12)
		cout << " AM" << endl;
	else 
		cout << " PM" << endl;
}


int Time::hour() const
{
	return dm_hour;
}

int Time::minute() const
{
	return dm_minute;
}
