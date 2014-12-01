#include <iostream>			// so we can use 'strcpy()'
using namespace std;	

// so compiler knows class declaration (what methods and data members to expect)
#include "house_pet.h"	


/**
The method code (definition) can be done any where in the file, usually done
in a seperate .cpp file. The class decleration (its name, methods and data fields)
are usually placed in it's own .h header file and included when needed.
**/

// method names must be scoped to the class they belong

// default constructor
HousePet::HousePet()
{
	strcpy(dm_animalType,"unkown animal type");
	strcpy(dm_name, "no name");
	dm_age = 0;
}

// copy constructor
HousePet::HousePet( const HousePet &pet )
{
	strcpy(dm_animalType, pet.dm_animalType);
	strcpy(dm_name, pet.dm_name);
	dm_age = pet.dm_age;
}

// overloaded constructor
HousePet::HousePet( char *name )
{
	strcpy(dm_animalType,"unkown animal type");
	strcpy(dm_name, name);
	dm_age = 0;
}

// uses default arguments
// note, that default values are not included in method definition
// only in method declaration
HousePet::HousePet( char *type, char *name, int age)
{
	strcpy(dm_animalType,type);
	strcpy(dm_name, name);
	dm_age = age;
}

HousePet::~HousePet()
{
	// does nothing for now, but defined just for kicks
}


void HousePet::set_animal_type(char *type)
{
	strcpy(dm_animalType,type);
}

void HousePet::set_name(char *name)
{
	strcpy(dm_name,name);
}

