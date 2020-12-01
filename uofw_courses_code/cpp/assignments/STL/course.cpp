/***
*	Method Code for the Course class.  
***/

#include <iostream>
using namespace std;

#include "course.h"
#include "student.h"

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
	dm_name = "Unknown Course Name";
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
	dm_studentList = c.dm_studentList;

	return *this;
}

// To Print out the data members (Course names and list of students in each).
void Course::print() // const
{
	Student *s;	// Local pointer to student object.
	int i = 1;

	cout << "Course Name: " << dm_name << endl;
	
	if (dm_studentList.empty())
	{
		cout << "Course has no Students\n" << endl;
		return;
	}

	// Display list of Students in the Course.
	cout << "The following Students are enrolled in this Course:" << endl;
	
	list<Student*>::iterator itr = dm_studentList.begin();

	while (itr != dm_studentList.end())
	{
		s = *itr;  // Get a Student from the list.

		cout << i << ". - " << (*s);
		cout << endl;
		i++;
		itr++;
	}
	cout << endl;
}

// To print out the list of students in a Course only.
void 
Course::print_students () // const
{
	list<Student*>::iterator itr = dm_studentList.begin();
	Student *s;  // Local pointer to a student object.
	int i = 1;

	if (dm_studentList.empty())
	{
		cout << "\nCourse has no Students\n" << endl;
		return;
	}

	cout << "\nThe following students are enrolled in that course:" << endl;
	
	

	while (itr != dm_studentList.end())
	{
		s = *itr;  // Get a student from the list.

		cout << i << ". - " << (*s) << " - Student #" << s->get_student_number();
		cout << endl;
		i++;
		itr++;
	}
}



// To maintain a list of all students enrolled in a particular course.
bool Course::add_student(Student *s)
{
	list<Student*>::iterator itr = dm_studentList.begin();

	// Check for bad data.
	if (s == 0)  
		return false;

	// Check to see if that student is already in the course's list.
	if (!dm_studentList.empty())  // no point if the list is empty.
	{
		// Use Object Address Comparison to see if the selected student is 
		// already in that course.
		while (itr != dm_studentList.end())
		{
			if (*itr == s)
			{
				cerr << "\n*** Course already contains that student ***" << endl;
				return false;
			}
			itr++;
		}
	}

	dm_studentList.push_back(s);  // Add the course to the list.
	return true;
}


bool Course::remove_student(Student *s)
{
	list<Student*>::iterator itr = dm_studentList.begin();
	list<Student*>::iterator next;

	if (s == 0)  // Check for bad data.
		return false;

	// Check if the student to be removed is in the list.
	if (!dm_studentList.empty())  // Make sure the list is not empty.
	{

		while (itr != dm_studentList.end())
		{
			// Use Object Address Comparison to remove the target course.
			if (*itr == s)
			{
				// Remove course from list and Update the linked list 
				next = dm_studentList.erase(itr);

				// Do not delete the memory associated with this pointer now because
				// the course does not OWN the student objects.  It has a list of 
				// nodes which point to students in the master student list.  Therefore
				// numerous other course could be using those same student objects.
				
				// Exit function once the target student has been found and removed.
				return true;
			}
			itr++;
		}
	}
	// If the course was not found
	return false;
}