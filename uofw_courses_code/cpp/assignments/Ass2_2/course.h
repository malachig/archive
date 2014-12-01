#ifndef _COURSE_H_
#define _COURSE_H_

#include <iostream>
using namespace std;

#include "my_string.h"


class Course; // forward class declaration

// Apparently microsoft C++ wants this function declared outside of the friend
// class.  Since one of it's arguments is a 'Course' reference, we must also have 
// the forward class declaration above.

ostream& operator<<(ostream &os, const Course &c);

/***
*	Simple class that keeps track of course information.  For now this 
*   is simply a string object containing the course name.
*   Has an overloaded constructor that accepts the course name.  Has method
*   to print out the course information.  The assignment operator '=' is 
*   overloaded.  The '<<' operator is also overloaded to print the course 
*	information.  Accessor methods: set_course() allows user to set the course name
*   and get_course() returns a constant reference to the course name.
***/

class Course
{
	// Overloading the '<<' operator
	friend ostream& operator<<(ostream &os, const Course &c);

public:
	Course();					// Default constructor.
	Course(const Course &c);	// Copy constructor.
	Course(const String &name); // Overloaded constructor.

	~Course();	// Destructor

	const Course& operator=(const Course &c);

	// Accessors
	const String & get_name() const {return dm_name;};
	void set_name(const String &name) {dm_name.set_data(name);}

	void print() const;


private:
	// DATA MEMBERS 
	String dm_name;

};

#endif
