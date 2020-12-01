/*
*  Date: Jan. 13 2002 
*  Author: Malachi Griffith
*  Purpose: This program introduces the concept of a simple class and 
*  a main() program with utilizes that class.
*/

#include <iostream>
using namespace std;

/***
	class: HousePet

	This class represents a pet that can be kept in the house.
	Its attributes are type, name and age.
***/

class HousePet
{
// The *public* interface, the methods (and sometimes data) that external 
// users of this class can access at any time.
public:
	// mutators, change the state of an object. Function Prototypes:
	void set_animal_type(char *type);
	void set_name(char *name);
	void set_age(int age);

	// accessors, allow access to object properties. Function Prototypes:
	char* get_animal_type();
	char* get_name();
	int get_age();

	// data members should be *private*, ENCAPSULATION & DATA HIDING!!
private:
	char dm_animalType[100];
	char dm_name[100];
	int dm_age;
};	// <=== note the semi-colon ';' at the end of the class block.

/*
	main() function *uses* the class declared above.
*/
int main()
{
	HousePet theCat;	// create an OBJECT
	HousePet theDog;	// create another object

	// set object state information, i.e. the data fields.
	theCat.set_animal_type("Felis Horribilis");
	theCat.set_name("Garfield");
	theCat.set_age(11);

	theDog.set_animal_type("Canis gigantis");
	theDog.set_name("Snoopy");
	theDog.set_age(22);

	// access the object's state information via accessor methods.
	cout << theCat.get_name() << " is ";
	cout << theCat.get_age() << " years old." << endl;
	cout << theDog.get_name() << " is ";
	cout << theDog.get_age() << " years old." << endl;

	return 0;
}

/**
The method code (definition) can be done anywhere in the file, usually done 
in a seperate .cpp file.  The class declaration (its name, methods and data
fields) are usually placed in it's own .h header file and included where
needed.
**/

// method names must be scoped to the class they belong to;

void HousePet::set_animal_type(char *type)
{
	strcpy(dm_animalType, type);
}

void HousePet::set_name(char *name)
{
	strcpy(dm_name, name);
}

void HousePet::set_age(int age)
{
	dm_age = age;
}

char* HousePet::get_animal_type()
{
	return dm_animalType;
}

char* HousePet::get_name()
{
	return dm_name;
}

int HousePet::get_age()
{
	return dm_age;
}



