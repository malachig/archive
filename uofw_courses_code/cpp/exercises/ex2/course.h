#ifndef _COURSE_H_  // To avoid multiple and recursive inclusions

#define _COURSE_H_

/***
 * This file defines the class Course.  -> Class declaration.
 * The "Course" class encapsulates the following information:
 * 	Course number <int>
 * 	Course name <char array> of max size 100.
 * The required mutators and accessors are provided.
***/

class Course 
{
	// The *public interface, the methods (and sometimes data) that 
	// external users of this class can access at any time.

public:
	// Constructors and destructors
	Course();  // default constructor.

	Course(const Course &crs);
	Course(char *name); // overloaded constructor.
	Course(char *name, double number = 0); // uses default arguments.

	// Mutators, change the state of an object.
	void set_course_name(char *name);
	void set_course_number(double number); 

	// Accessors, allow access to object properties
	char* get_course_name() {return dm_courseName;}
	double get_course_number() {return dm_number;}

	// Data members should be *private, ENCAPSULATION & DATA HIDING!!

private:
	char dm_courseName[100];
	double dm_number;

};  // <=== Note the semi-colon ';' at the end of the class block.

#endif
	





















