/***
*	The method code (definition) for the Person Class.
*	Located in a seperate .cpp file for clarity and modularity. 
*	The class decleration (its name, methods and data fields) are found
*	in the person.h header file and included when needed.
***/

// So compiler knows class declaration (methods and data members to expect)
#include "person.h"	

// Note: Method names must be scoped to the class they belong.

// Overloadeding the '<<' Operator
ostream& operator<<(ostream &os, const Person &p)
{
	os << p.dm_name;
	return os;
}


// Default Constructor.
Person::Person()
{
	dm_name = "Unknown Name";
	dm_addr = "Unknown Address";
}


// Copy Constructor, uses data member initialization list
Person::Person(const Person &p)
: dm_name(p.dm_name), dm_addr(p.dm_addr)
{
  
}


// Overloaded Constructor.
// This one allows user to initialze Person attributes when the Person
// object is created.
Person::Person(const String &name, const String &addr)
: dm_name(name), dm_addr(addr)
{
  
}


Person::~Person()
{
	// Does nothing for now, but defined just for kicks.
}


// Overloading the assignment '=' operator.
const Person& Person::operator=(const Person &p)
{
	if (this == &p)		// If LHS is the same as RHS
		return *this;	
	
	dm_name = p.dm_name;
	dm_addr = p.dm_addr;
	
	return *this;
}


// Method Code for Mutators

void 
Person::set_name(const String &name)
{
	dm_name = name;
}


void 
Person::set_address(const String &addr)
{
	dm_addr = addr;
}


// Method Code for Accessors

void 
Person::print() const
{
	// Now actually display the info to the screen.
	cout << "\nPERSONAL INFO" << endl;
	cout << "Name: " << dm_name << endl;
	cout << "Address: " << dm_addr << endl;
	cout << endl;
}

