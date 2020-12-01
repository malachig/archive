/*
*  Date: Jan. 13 2002 
*  Author: Malachi Griffith
*  Purpose: This program introduces the concept of a simple class and 
*  a main() program with utilizes that class.
*/

#include <iostream>
#include <string>
using namespace std;

/***
	class: Student 

	This class represents a university student.
	Its attributes are name and student number.
***/

class Student 
{
// The *public* interface, the methods (and sometimes data) that external 
// users of this class can access at any time.
public:
	// mutators, change the state of an object. Function Prototypes:
	void set_name();
	void set_number();

	// accessors, allow access to object properties. Function Prototypes:
	string get_name();
	int get_number();

	// data members should be *private*, ENCAPSULATION & DATA HIDING!!
private:
	string dm_name;
	int dm_number;
};	// <=== note the semi-colon ';' at the end of the class block.

/***
	class: Course

	This class represents a university course.
	Its attributes are name and course number.
***/

class Course
{
// The *public* interface, the methods (and sometimes data that external
// users of this class can access at any time.
public:
	// mutators, change the state of an oblect. Function Prototypes:
	void set_cr_name(char *cr_name);
	void set_cr_number(float cr_number);

	// accessors, allow access to oblect propeties. Function Prototypes:
	char* get_cr_name();
	float get_cr_number();

	// data members should be *private*, ENCAPSULATION & DATA HIDING!!
private:
	char dm_cr_name[100];
	float dm_cr_number;
};


/*
	main() function *uses* the class declared above.
*/
int main()
{
	Student student1;	// create an OBJECT
	Student student2;	// create another object
	Course course1;
	Course course2;

	// set object state information, i.e. the data fields.
	student1.set_name();
	student1.set_number();
	student2.set_name();
	student2.set_number();
	
	course1.set_cr_name("Intro to C++");
	course1.set_cr_number(91.2947);	
	course2.set_cr_name("Programming in C");
	course2.set_cr_number(91.1905);

	// access the object's state information via accessor methods.
	cout << student1.get_name() << " has ";
	cout << student1.get_number() << " for a student number." << endl;
	cout << student2.get_name() << " has ";
	cout << student2.get_number() << " for a student number." << endl;

	cout << course1.get_cr_name() << " has ";
	cout << course1.get_cr_number() << " as its number." << endl;
	cout << course2.get_cr_name() << " has ";
	cout << course2.get_cr_number() << " as its number." << endl;

	return 0;
}

/**
The method code (definition) can be done anywhere in the file, usually done 
in a seperate .cpp file.  The class declaration (its name, methods and data
fields) are usually placed in it's own .h header file and included where
needed.
**/

// method names must be scoped to the class they belong to;

void Student::set_name()
{
	string name_f, name_l;
	
	cout << "Please enter the first name: ";
	cin >> name_f;
	cout << "Please enter the last name: ";
	cin >> name_l;
	dm_name = name_f + " " + name_l;
}

void Student::set_number()
{
	int number;

	cout << "Please enter the student number as an integer: ";
	cin >> number;
	dm_number = number;
}

string Student::get_name()
{
	return dm_name;
}

int Student::get_number()
{
	return dm_number;
}


// Now for the Course Class methods:

void Course::set_cr_name(char *cr_name)
{
        strcpy(dm_cr_name, cr_name);
}
 
void Course::set_cr_number(float cr_number)
{
        dm_cr_number = cr_number;
}
 
char* Course::get_cr_name()
{
        return dm_cr_name;
}
 
float Course::get_cr_number()
{
        return dm_cr_number;
}

