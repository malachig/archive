#ifndef _STUDENT_H_ // To avoid multiple and recursive inclusions
#define _STUDENT_H_

#include "my_string.h"   // Student name and address will be string objects.
#include "linked_list.h" // To keep a linked list of courses for each student.
#include "course.h"		 // Describes the course objects in the afore mentioned list.

#include <iostream>	 // Needed for print() behaviour.
using namespace std;

/***
*	class: Student
*
*	This class represents a student at a university.  Students have the 
*	following basic attributes:
*		student number (unsigned long).
*		student name (String, class provided by Rodrigo Vivanco).
*		address (String as well).
*		course list (Course object).
*	Behaviours:
*		Various constructors.
*		Set name, address, or number.
*		Get name, address, or number.
*	    Add course.
*		Remove course.
*		Get list of courses.
*		Print student information.
***/

class Student
{
// The *public* interface, the methods (and sometimes data) that external
// users of this class can access at any time.
public:

  // Constructors and Destructors.
	
  Student();                  // Default constructor.
  Student(const Student &s);  // Copy constructor.
  
  // Overloaded constructor.
  Student(const String &name, const String &addr, unsigned long number);

  ~Student();

  // Mutators, change the state of an object.
  void  set_student_name(const String &name);
  void  set_student_address (const String &addr);
  void  set_student_number(unsigned long number) {dm_studentNumber = number;}

  // Accessors, allow access to object properties.
  // Do not allow name or address to be changed and return as a reference.
  const String & get_student_name() const {return dm_studentName;}
  const String & get_student_addr() const {return dm_studentAddr;}
  unsigned long get_student_number() const;

  // Print accessor.  Prints student object to screen.
  void print () const;
  
  // Print student's courses only.
  void print_courses () const;

  // Methods to deal with the student's courses.
  LinkedList & get_course_list() {return dm_courseList;}
  void add_course(Course *c);
  void remove_course(Course *c);
  
// Data members should be *private*, ENCAPSULATION & DATA HIDING!!
private:
	
  String dm_studentName;
  String dm_studentAddr;
  unsigned long dm_studentNumber;
  LinkedList dm_courseList;  // Linked list of courses for the student object.
};	

#endif