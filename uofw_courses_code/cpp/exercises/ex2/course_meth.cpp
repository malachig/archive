/***
 * This is the method code definition for the class Course.  In other words
 * each constructor and function found in the Course class is actually 
 * defined in here.
***/

#include <iostream>
using namespace std;

// So compiler knows class declaration (what methods and data members are 
// to be expected:
#include "course.h"

/***
 * The method code follows
***/

// Default constructor
Course::Course()
{
	strcpy(dm_courseName, "Unknown Course");
	dm_number = 0;
}

// Copy constructor - ie. if a course is already defined, copy its values
// to a new course by sending the course variable name to the constructor.
Course::Course(const Course &crs)
{
	strcpy(dm_courseName, crs.dm_courseName);
	dm_number = crs.dm_number;
}

// Overloaded constructor
Course::Course(char *name)
{
	strcpy(dm_courseName, name);
	dm_number = 0;
}

// uses default arguments.
// note, that default values are not included in method definition
// only in method declaration.
Course::Course(char *name, double number)
{
	strcpy(dm_courseName, name);
	dm_number = number;
}

// Now begin the actual functions 
void Course::set_course_name(char *name)
{
	strcpy(dm_courseName, name);
}

void Course::set_course_number(double number)
{
	dm_number = number;
}








