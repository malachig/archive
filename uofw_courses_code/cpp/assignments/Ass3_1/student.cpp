/***
*	The method code (definition) for the Student Class.
*	Located in a seperate .cpp file for clarity and modularity. 
*	The class decleration (its name, methods and data fields) are found
*	in the student.h header file and included when needed.
***/

// So compiler knows class declaration (methods and data members to expect)
#include "student.h"	

// Note: Method names must be scoped to the class they belong.

// Overloadeding the '<<' Operator
ostream& operator<<(ostream &os, const Student &s)
{
	os << s.dm_studentName;
	return os;
}


// Default Constructor.
Student::Student()
{
  dm_studentName.set_data("Unknown Name");
  dm_studentAddr.set_data("Unknown Address");
  dm_studentNumber = 0;
}


// Copy Constructor, uses data member initialization list
Student::Student(const Student &s)
: dm_studentName(s.dm_studentName), dm_studentAddr(s.dm_studentAddr)
{
  dm_studentNumber = s.dm_studentNumber;
}


// Overloaded Constructor.
// This one allows user to initialze Student attributes when the student
// object is created.
Student::Student(const String &name, const String &addr, unsigned long number)
: dm_studentName(name), dm_studentAddr(addr)
{
  dm_studentNumber = number;
}


Student::~Student()
{
  // Does nothing for now, but defined just for kicks.
}


// Overloading the assignment '=' operator.
const Student& Student::operator=(const Student &s)
{
	if (this == &s)		// If LHS is the same as RHS
		return *this;	
	
	dm_studentName = s.dm_studentName;
	dm_studentAddr = s.dm_studentAddr;
	dm_studentNumber = s.dm_studentNumber;
	dm_courseList = s.dm_courseList;

	return *this;
}


// Method Code for Mutators

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


// Method Code for Accessors

unsigned long 
Student::get_student_number() const
{
	return dm_studentNumber;
}


void 
Student::print() const
{
	Course *c;				// Local pointer to course object.

	// Now actually display the info to the screen.
	cout << "\nThe student's name is: " << dm_studentName << endl;
	cout << "The student's address is: " << dm_studentAddr << endl;
	cout << "The student # is: " << dm_studentNumber << endl;

	if (dm_courseList.size() < 1)
	{
		cout << "\nStudent has no courses\n" << endl;
		return;
	}

	// Display list of courses for the student.
	cout << "\nThe student is enrolled in the following courses:" << endl;
	
	for (int i = 0; i < dm_courseList.size(); i++)
	{
		c = (Course*)dm_courseList[i];  // Get a course from the list.

		cout << i+1 << ". - " << (*c);
		cout << endl;
	}
	cout << endl;
}


void 
Student::print_courses () const
{
	Course *c;  // Local pointer to a course object.

	cout << "That student is enrolled in the following courses:" << endl;
	
	for (int i = 0; i < dm_courseList.size(); i++)
	{
		c = (Course*)dm_courseList[i];  // Get a course from the list.

		cout << i+1 << ". - " << (*c);
		cout << endl;
	}
}


void 
Student::add_course(Course *c)
{
	// Check for bad data.
	if (c == 0)  
		return;

	// Check to see if that course is already in the student's list.
	if (dm_courseList.size() > 0)  // no point if the list is empty.
	{
		// Use Object Address Comparison to see if the selected course is 
		// already in the student's list.
		for (int i = 0; i < dm_courseList.size(); i++)
		{
			if (dm_courseList[i] == c)
			{
				cerr << "\n*** Student already enrolled in that course ***" << endl;
				return;
			}
		}
	}

	dm_courseList.add(c);  // Add the course to the list.
	cout << "Course Added" << endl;
}


void 
Student::remove_course(Course *c)
{
	Course *removeCourse;
	
	if (c == 0)  // Check for bad data.
		return;

	// Check if the course to be removed is in the list.
	if (dm_courseList.size() > 0)  // Make sure the list is not empty.
	{
		for (int i = 0; i < dm_courseList.size(); i++)
		{
			// Use Object Address Comparison to remove the target course.
			if (dm_courseList[i] == c)
			{
				// Remove course from list and Update the linked list 
				removeCourse = (Course*)dm_courseList.remove(i);

				// Do not delete the memory associated with this pointer now because
				// the student does not OWN the course objects.  It has a list of 
				// nodes which point to courses in the master course list.  Therefore
				// numerous other students could be using those same course objects.
				
				// Exit function once the target course has been found and removed.
				return;
			}
		}
		// If the course was not found display an error message
		cerr << "\n*** Student not enrolled in that course ***" << endl;
	}
}
