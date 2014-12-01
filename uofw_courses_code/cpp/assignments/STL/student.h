#ifndef _STUDENT_H_ // To avoid multiple and recursive inclusions
#define _STUDENT_H_


#include <iostream>	 // Needed for print() behaviour.
#include <list>
using namespace std;

#include "course.h"		 // Describes the course objects in the afore mentioned list.
#include "person.h"      // Needed to inherit attributes of base class person.

/***
*	class: Student - Inherits from the base class Person
*
*	This class represents a student at a university.  Students have the 
*	following basic attributes:
*		student number (unsigned long).
*		course list (Course object).
*	Behaviours:
*		Various constructors.
*		Set student number.
*		Get student number.
*	    Add course.
*		Remove course.
*		Get list of courses.
*		Print student information.
***/


class Student : public Person
{

// The *public* interface, the methods (and sometimes data) that external
// users of this class can access at any time.
public:

	// Constructors and Destructors.
	Student();                  // Default constructor.
	Student(const Student &s);  // Copy constructor.
  
	// Overloaded constructor.
	Student(unsigned long number);

	virtual ~Student();

	// Overloading the assignment operator
	const Student& operator=(const Student &s);

	// Mutators, change the state of an object.
	void  set_student_number(unsigned long number) {dm_studentNumber = number;}

	// Accessors, allow access to object properties.
	unsigned long get_student_number() const;

	// Print accessor.  Prints student object to screen.
	virtual void print (); // const;
  
	// Print student's courses only.
	void print_courses (); // const; 

	// Methods to deal with the student's courses.
	list<Course*>& get_course_list() {return dm_courseList;}
	void add_course(Course *c);
	void remove_course(Course *c);
  
// Data members should be *protected*, Allow subclasses to access.
protected:
	
	unsigned long dm_studentNumber;
	list<Course*> dm_courseList;  // STL list of courses for the student object.
};	


/***
*	class: GradStudent - A specialized class of Student
*
*	This class represents a Grad student.  Grad Students have the 
*	following basic attributes:
*		Thesis title
*		Advisors name.
*	Behaviours:
*		Various constructors.
*		Set thesis title or advisors name.
*		Get thesis title or advisors name.
*		Print Grad student information.
***/

class GradStudent : public Student
{
public:

	// Constructors and Destructors.
	GradStudent();                  // Default constructor.
	GradStudent(const GradStudent &gs);  // Copy constructor.
  
	// Overloaded constructor.
	GradStudent(const String &thesis, const String &advisor);
	
	virtual ~GradStudent();

	// Overloading the assignment operator
	const GradStudent& operator=(const GradStudent &gs);

	// Mutators, change the state of an object.
	void set_thesis_title(const String &thesis);
	void set_advisor_name(const String &advisor);
	
	// Accessors, allow access to object properties.  Do not allow thesis
	// title or advisor name to be changed by get method.
	const String & get_thesis_title() const {return dm_thesisTitle;}
	const String & get_advisor_name() const {return dm_advisorName;}

	// Print accessor.  Prints student object to screen.
	virtual void print (); // const;
  
// Data members should be *protected*, Allow subclasses to access.
protected:
	String dm_thesisTitle;
	String dm_advisorName;
};




/***
*	class: UnderGradStudent - A specialized class of Student
*
*	This class represents an UnderGrad student.  UnderGrads have the 
*	following basic attributes:
*		Project title
*	Behaviours:
*		Various constructors.
*		Set project title.
*		Get project title.
*		Print UnderGrad student information.
***/

class UnderGradStudent : public Student
{
public:

	// Constructors and Destructors.
	UnderGradStudent();                  // Default constructor.
	UnderGradStudent(const UnderGradStudent &ugs);  // Copy constructor.
  
	// Overloaded constructor.
	UnderGradStudent(const String &project);
	
	virtual ~UnderGradStudent();

	// Overloading the assignment operator
	const UnderGradStudent& operator=(const UnderGradStudent &ugs);

	// Mutators, change the state of an object.
	void set_project_title(const String &project);
	
	// Accessors, allow access to object properties.  Do not allow project
	// title to be changed by get method.
	const String & get_project_title() const {return dm_projectTitle;}

	// Print accessor.  Prints student object to screen.
	virtual void print (); // const;
  
// Data members should be *protected*, Allow subclasses to access.
protected:
	String dm_projectTitle;
};

#endif