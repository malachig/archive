#ifndef _PERSON_H_ // To avoid multiple and recursive inclusions
#define _PERSON_H_

#include "my_string.h" // Student name and address will be string objects.

#include <iostream>	 // Needed for print() behaviour.
using namespace std;

/***
*	class: Person
*
*	This class represents a Person.  People have the 
*	following basic attributes:
*		name (String, class provided by Rodrigo Vivanco).
*		address (String as well).
*	Behaviours:
*		Various constructors.
*		Set name or address.
*		Get name or address.
*		Print person information.
***/

class Person; // forward class declaration

// Apparently microsoft C++ wants this function declared outside of the friend
// class.  Since one of it's arguments is a 'Student' reference, we must also have 
// the forward class declaration above.

// Overloading the '<<' operator
ostream& operator<<(ostream &os, const Person &p);

class Person
{
	// Making this function a friend so that it can access data members directly
	friend ostream& operator<<(ostream &os, const Person &p);

public:
	
	// Constructors and Destructors.
	Person();                  // Default constructor.
	Person(const Person &p);  // Copy constructor.
  
	// Overloaded constructor.
	Person(const String &name, const String &addr);

	virtual ~Person();

	// Overloading the assignment operator
	const Person& operator=(const Person &p);

	// Mutators, change the state of an object.
	void  set_name(const String &name);
	void  set_address (const String &addr);
	
	// Accessors, allow access to object properties.
	// Do not allow name or address to be changed and return as a reference.
	const String & get_name() const {return dm_name;}
	const String & get_addr() const {return dm_addr;}
	
	// Print accessor.  Prints student object to screen.
	virtual void print () const;
  
	  
// Data members should be *protected*, Allow subclasses to access.
protected:
	
	String dm_name;
	String dm_addr;
};

#endif