#include <iostream>
using namespace std;

int main()
{
	// declare variables
	short	shoeSize;
	int		age;

	char	firstName[100];	// automaic array
	char	*lastName;			// dynamic array, must allocate storage

	cout << "Please enter the following information:" << endl;

	cout << "first name: ";
	cin >> firstName;

	lastName = new char[100];
	cout << "last name: ";
	cin >> lastName;

	cout << "age: ";
	cin >> age;

	cout << "shoe size: ";
	cin >> shoeSize;

	bool	canVote;	// can declare variables anywhere inside a block

	if ( age - shoeSize >= 11 )
		canVote = true;
	else
		canVote = false;

	cout << firstName << ' ' << lastName;
	cout << ": according to our calculations you ";

	// hey, clever use of '?', who'd have thunk it !!
	canVote ? cout << "can" : cout << "can't";

	cout << " vote." << endl;		// *** insert clever action here ***

	// return dynamic storage to the system heap
	delete [] lastName;

	return 0;
}