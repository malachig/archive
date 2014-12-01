#include <iostream>
using namespace std;

/***
	class: HousePet

	this class represents a pet that can be kept in a house, flies and other
	little insect creatures are likely not to be pets, but who knows what 
	weird tastes people have...
***/

class HousePet
{
// the *public* interface, the methods (and sometime data) that external
// users of this class can access at any time
public:

	// mutators, change the state of an object
	void	set_animal_type( char *type );
	void	set_name( char *name);
	void	set_age( int age );

	// accessors, allow access to object properties
	char*	get_animal_type();
	char*	get_name();
	int		get_age();


// data members should be *private*, ENCAPSULATION & DATA HIDING!!
private:
	
	char	dm_animalType[100];
	char	dm_name[100];
	int		dm_age;
};	// <==== note the semi-colon ';' at the end of the class block


/*
	main() function *uses* the class declared above.
*/
int main()
{
	HousePet	theCat;	// create an OBJECT
	HousePet	theDog;	// create another OBJECT

	// set object state information, i.e. the data fields
	theCat.set_animal_type("Felis");
	theCat.set_name("Garfield");
	theCat.set_age(11);

	theDog.set_animal_type("Canis");
	theDog.set_name("Snoopy");
	theDog.set_age(22);

	// access the object's state information via accessor methods
	cout << theCat.get_name() << " is ";
	cout << theCat.get_age() << " years old." << endl;
	cout << theDog.get_name() << " is ";
	cout << theDog.get_age() << " years old." << endl;

	return 0;
}



/**
The method code (definition) can be done any where in the file, usually done
in a seperate .cpp file. The class decleration (its name, methods and data fields)
are usually placed in it's own .h header file and included when needed.
**/

// method names must be scoped to the class they belong

void HousePet::set_animal_type(char *type)
{
	strcpy(dm_animalType,type);
}

void HousePet::set_name(char *name)
{
	strcpy(dm_name,name);
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