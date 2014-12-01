/***
*	The method code (definition) for the Student Class.
*	Located in a seperate .cpp file for clarity and modularity. 
*	The class decleration (its name, methods and data fields) are found
*	in the student.h header file and included when needed.
***/

// So compiler knows class declaration (methods and data members to expect)
#include "student.h"	

// Note: Method names must be scoped to the class they belong.

// Default Constructor.
Student::Student()
{
  dm_studentName.set_data("Unknown Name");
  dm_studentAddr.set_data("Unknown Address");
  dm_studentNumber = 1;
}

// Copy Constructor.
Student::Student(const Student &s)
{
  dm_studentName.set_data(s.dm_studentName);
  dm_studentAddr.set_data(s.dm_studentAddr);
  dm_studentNumber = s.dm_studentNumber;
}

// Overloaded Constructor.
// This one allows user to initialze Student attributes when the student
// object is created.
Student::Student(const String &name, const String &addr, unsigned long number)
{
  dm_studentName.set_data(name);
  dm_studentName.set_data(addr);
  dm_studentNumber = number;
}


Student::~Student()
{
  // Does nothing for now, but defined just for kicks.
}



// Method code for Mutators

void 
Student::set_student_name(const String &name)
{
	dm_studentName.set_data(name);
}


void 
Student::set_student_address(const String &addr)
{
	dm_studentAddr.set_data(addr);
}


unsigned long 
Student::get_student_number() const
{
	return dm_studentNumber;
}


void 
Student::print() const
{
	// Need some temporary character arrays because 'cout' does not like 
	// our string objects.
	const char *tempName;
	const char *tempAddr;

	// Get the values using string behaviour, get_data().
	tempName = dm_studentName.get_data();
	tempAddr = dm_studentAddr.get_data();

	// Now actually display the info to the screen.
	cout << "\nThe student's name is: " << tempName << endl;
	cout << "The student's address is: " << tempAddr << endl;
	cout << "The student # is: " << dm_studentNumber << endl;
}

