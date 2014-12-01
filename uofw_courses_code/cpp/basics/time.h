#ifndef _TIME_H_ // To avoid multiple and recursive inclusions
#define _TIME_H_

class Time 
{
public:
	Time();		// Constructor.

	// Public interface, member functions.
	void set_time(int hour, int minute=0, int second=0);
	void set_time(const Time &t);
	void print_military();
	void print_standard();

	// Accessor functions
	int hour() const;
	int minute() const;

private:
	// Data members.
	int dm_hour;	// 0-23
	int dm_minute;  // 0-59
	int dm_second;	// 0-59
};

/* Classes are define with the keyword class and delineated with { and }.
*  The class definition must be terminated with a semicolon.
*/

#endif