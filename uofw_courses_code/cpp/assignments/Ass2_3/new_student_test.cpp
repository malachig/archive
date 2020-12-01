/***
*	This program is a test function for the student class.
*	Several course objects will be created as a linked list within a student object.
*	Note: Normally I would not include so much code in main that need not be there,
*	however this is only a preliminary test to prove that course tracking features of 
*   my student class work.  ie. Can I ADD courses (checking to make sure that the course
*	has not already been added), REMOVE course (if the course is in the list), and 
*	PRINT the student object with the courses that student belongs to.
***/

#include <iostream>
using namespace std;

#include "student.h"

void main()
{
	Student s1;
	Student s2;

	Course *c1;
	Course *c2;
	Course *c3;

	s1.set_student_name("Malachi Griffith");
	s1.set_student_address("Winnipeg");
	s1.set_student_number(1056419);

	s2.set_student_name("Scott Ticknor");
	s2.set_student_address("Ottawa");
	s2.set_student_number(1111111);

	c1 = new Course("Biology");
	c2 = new Course("Chemistry");
	c3 = new Course("Anatomy");

	s1.add_course(c1);
	s1.add_course(c2);
	s1.add_course(c3);
	s1.add_course(c2); // try to add a course already in list.

	s2.add_course(c3);

	s1.print();
	s2.print();

	s1.remove_course(c2);
	s1.remove_course(c3);

	s1.print();
}

