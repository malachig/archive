#include <iostream>
using namespace std;  // so we can use 'strcpy()'

// so compiler knows class declaration (what methods and data members expected)
#include "house_pet.h"

/***
 * The method code follows
***/

// Default constructor
HousePet::HousePet()
{
	strcpy(dm_animalType, "Unknown animal type");
	strcpy(dm_name, "No name");
	dm_age = 0;
}

// Copy constructor - ie if a pet is already defined copy its values 
// to a new pet by sending the pet name to the constructor.
HousePet::HousePet(const HousePet &pet)
{
	strcpy(dm_animalType, pet.dm_animalType);
	strcpy(dm_name, pet.dm_name);
	dm_age = pet.dm_age;
}

// Overloaded constructor
HousePet::HousePet(char *name)
{
	strcpy(dm_animalType, "Unknown animal type");
	strcpy(dm_name, name);
	dm_age = 0;
}

// uses default arguments
// note, that default values are not included in method definition
// only in method declaration
HousePet::HousePet(char *type, char *name, int age)
{
	strcpy(dm_animalType, type);
	strcpy(dm_name, name);
	dm_age = age;
}

HousePet::~HousePet()
{
	// does nothing for now, but defined just for kicks
}

// Now begin the actual functions
void HousePet::set_animal_type(char *type)
{
	strcpy(dm_animalType, type);
}

void HousePet::set_name(char *name)
{
	strcpy(dm_name, name);
}

