#include <iostream>
using namespace std;

#include "course.h"

#define errMsg "<Error> Course::"

// Overloaded '<<' Operator
ostream& operator<<(ostream &os, const Course &c)
{
	os << c.dm_name;
	return os;
}

// Default Constructor.
Course::Course()
{
	dm_name.set_data("Unknown Course Name");
}

// Copy Constructor
Course::Course(const Course &c)
: dm_name(c.dm_name)
{

}

// Overloaded Constructor 
Course::Course(const String &name)
: dm_name(name)
{

}

// Destructor.  Included for future development
Course::~Course()
{

}

// Overloaded assignment '=' operator.
const Course& Course::operator=(const Course &c)
{
	if (this == &c)
		return *this;	// LHS is the same as RHS
	
	dm_name = c.dm_name;

	return *this;
}

// Print out the data members.
void Course::print() const
{
	cout << "Course Name: " << dm_name << endl;
}
