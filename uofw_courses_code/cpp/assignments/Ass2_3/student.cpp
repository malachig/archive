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
	Course *c;  // Local pointer to course object.

	// Get the values using string behaviour, get_data().
	tempName = dm_studentName.data();
	tempAddr = dm_studentAddr.data();

	// Now actually display the info to the screen.
	cout << "\nThe student's name is: " << tempName << endl;
	cout << "The student's address is: " << tempAddr << endl;
	cout << "The student # is: " << dm_studentNumber << endl;

	cout << "The student is enrolled in the following courses:" << endl;
	
	for (int i = 0; i < dm_courseList.size(); i++)
	{
		c = (Course*)dm_courseList[i];  // Get a course from the list.

		cout << i << ". - " << (*c);
		cout << endl;
	}
}

void 
Student::add_course(Course *c)
{
	Course *tempCourse;
	int test;

	if (c == 0)  // Check for bad data.
		return;


	// Check to see if that course is already in the list.
	if (dm_courseList.size() > 0)  // no point if the list is empty.
	{
		for (int i = 0; i < dm_courseList.size(); i++)
		{
			tempCourse = (Course*)dm_courseList[i];
			
			test = strcmp((tempCourse->get_name()).data(), (c->get_name()).data());

			if (test == 0)
			{
				cerr << "Student already enrolled in that course!" << endl;
				return;
			}
		}
		
	}
	dm_courseList.add(c);  // Add the course to the list.
}



void 
Student::remove_course(Course *c)
{
	Course *tempCourse;
	Course *removeCourse;
	int test;

	if (c == 0)  // Check for bad data.
		return;

	// Check if the course to be removed is in the list.
	if (dm_courseList.size() > 0)  // Can't be in list if it is empty.
	{
		for (int i = 0; i < dm_courseList.size(); i++)
		{
			tempCourse = (Course*)dm_courseList[i];

			test = strcmp((tempCourse->get_name()).data(), (c->get_name()).data());

			if (test == 0)
			{
				// Remove course from list and Update the linked list 
				removeCourse = (Course*)dm_courseList.remove(i);
				
				// Free memory used by that course.
				delete removeCourse;
				
				// Exit function once the target course has been found and removed.
				return;
			}
		}
	// ??? Include error message if the course is not found ???
	}
}
