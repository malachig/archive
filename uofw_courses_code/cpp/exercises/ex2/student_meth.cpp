/***
 * This is the method code definition for the class Student.  In other words
 * each constructor and function found in the Student class is actually 
 * defined in here.
***/

#include <iostream>
#include <string>
using namespace std;

// So compiler knows class declaration (what methods and data members are 
// to be expected:
#include "student.h"

/***
 * The method code follows
***/

// Default constructor
Student::Student()
{
	strcpy(dm_studentName, "Unknown Student");
	dm_number = 0;
}

// Copy constructor - ie. it a student is already defined, copy its values
// to a new student by sending the student variable name to the constructor.
Student::Student(const Student &stu)
{
	strcpy(dm_studentName, stu.dm_studentName);
	dm_number = stu.dm_number;
}

// Overloaded constructor
Student::Student(char *name)
{
	strcpy(dm_studentName, name);
	dm_number = 0;
}

// uses default arguments.
// note, that default values are not included in method definition
// only in method declaration.
Student::Student(char *name, int number)
{
	strcpy(dm_studentName, name);
	dm_number = number;
}

// Now begin the actual functions 
void Student::set_student_name(char *name)
{
	strcpy(dm_studentName, name);
}

void Student::set_student_number(int number)
{
	dm_number = number;
}






