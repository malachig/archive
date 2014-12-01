#ifndef _HOUSE_PET_H_		// to avoid multiple and recursive inclusions
#define _HOUSE_PET_H_


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

	// constructors and destructors
	
	HousePet();												// default constructor
	HousePet( const HousePet &pet );	// copy constructor
	HousePet( char *name );						// overloaded constructor
	HousePet( char *type, char *name, int age=1);	// uses default arguments

	~HousePet();

	// mutators, change the state of an object
	void	set_animal_type( char *type );
	void	set_name( char *name);
	void	set_age( int age ) { dm_age = age; }	// inline method

	// accessors, allow access to object properties
	char*	get_animal_type() { return dm_animalType; }
	char*	get_name() { return dm_name; }
	int		get_age() { return dm_age; }


// data members should be *private*, ENCAPSULATION & DATA HIDING!!
private:
	
	char	dm_animalType[100];
	char	dm_name[100];
	int		dm_age;
};	// <==== note the semi-colon ';' at the end of the class block


#endif