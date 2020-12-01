/***
*	Method Code for the Course class.  
***/

#include <iostream>
using namespace std;

#include "course.h"

// Include this incase I want to use it later.
#define errMsg "<Error> Course::"

// Overloadeding the '<<' Operator
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

// Overloading the assignment '=' operator.
const Course& Course::operator=(const Course &c)
{
	if (this == &c)		// If LHS is the same as RHS
		return *this;	
	
	dm_name = c.dm_name;

	return *this;
}

// To Print out the data members (Course names for now).
void Course::print() const
{
	cout << "Course Name: " << dm_name << endl;
}
