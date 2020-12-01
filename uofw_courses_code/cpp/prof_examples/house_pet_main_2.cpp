#include <iostream>
using namespace std;

#include "house_pet.h"	// unless using namespaces, use 'tradition' #include

// global constants
const int MAX_NUM_PETS = 100; // at most 100 pets, should be plenty for now


// global variables, try to minimize these
HousePet	thePetsG[MAX_NUM_PETS];	// master pet list
int				numPetsG = 0;						// total number of pets in the master list


// menu choices available to user
enum MenuChoiceEnum { ADD_PET, REMOVE_PET, PRINT_PET_LIST, QUIT };


// *** function prototypes ***

// prints menu options, returns user input
MenuChoiceEnum get_menu_choice();

// carries out menu action entered by user
void perform_menu_action( MenuChoiceEnum menuChoice );
void add_pet();
void remove_pet();
void print_pet_list();






int main()
{
	MenuChoiceEnum	menuChoice;

	do
	{
		menuChoice = get_menu_choice();
		perform_menu_action( menuChoice );
	}
	while (menuChoice != QUIT);

	return 0;
}


// prints menu options, returns user input
MenuChoiceEnum get_menu_choice()
{
	MenuChoiceEnum menuEntryNum = QUIT;
	char	ch;	// if use int type for input, cin chokes on a letter	

	cout << "\n\n\n\n";	
	do
	{
		cout << "Enter choice (0 - 3 ): \n";
		cout << "0 - ADD PET\n"; 
		cout << "1 - REMOVE PET PET\n"; 
		cout << "2 - PRINT PET LIST\n"; 
		cout << "3 - QUIT\n"; 

		cin >> ch;			// can't put into user defined type (yet)
		ch = ch - '0';	// convet ascii 'number' into digital number
		menuEntryNum = (MenuChoiceEnum)ch; // now assign it to the user defined type
	} 
	while ( menuEntryNum<ADD_PET || menuEntryNum>QUIT );	// repeat until valid entry

	return menuEntryNum;
}

// as the name suggests, carry out menu action
void perform_menu_action( MenuChoiceEnum menuChoice )
{
	switch( menuChoice )
	{
	case ADD_PET:
		add_pet();
		break;
	
	case REMOVE_PET:
		remove_pet();
		break;

	case PRINT_PET_LIST:
		print_pet_list();
		break;
	}

}

// ** the following functions manipulate the master pet list ** //

void add_pet()
{
cerr << "*** add_pet() *** "<<endl;
}

void remove_pet()
{
cerr << "*** remove_pet() *** "<<endl;
}

void print_pet_list()
{
cerr << "*** print_pet_list() *** "<<endl;
}