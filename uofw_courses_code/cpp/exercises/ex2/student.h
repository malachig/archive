#ifndef _STUDENT_H_  // To avoid multiple and recursive inclusions

#define _STUDENT_H_

/***
 * This file defines the class Student.  -> Class declaration.
 * The "Student" class encapsulates the following information:
 * 	Student number <int>
 * 	Student name <char array> of max size 100.
 * The required mutators and accessors are provided.
***/

class Student
{
	// The *public interface, the methods (and sometimes data) that 
	// external users of this class can access at any time.

public:
	// Constructors and destructors
	Student();  // default constructor.

	Student(const Student &stu);
	Student(char *name); // overloaded constructor.
	Student(char *name, int number = 0); // uses default arguments.

	// Mutators, change the state of an object.
	void set_student_name(char *name);
	void set_student_number(int number); 

	// Accessors, allow access to object properties
	string get_student_name() {return dm_studentName;}
	int get_student_number() {return dm_number;}

	// Data members should be *private, ENCAPSULATION & DATA HIDING!!

private:
	char dm_studentName[100];
	int dm_number;

};  // <=== Note the semi-colon ';' at the end of the class block.

#endif
	





















