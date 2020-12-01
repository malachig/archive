/*
*  Date: January 13 2002 
*  Author: Malachi Griffith
*  Purpose: Illustrate the basic usage of input and coutput functions:
*  cout and cin, as well as the creation and deletion of some dynamic 
*  memory.
*/

#include <iostream>
using namespace std;

int main()
{
	// declare variables
	short shoeSize;
	int age;

	char firstName[100];	// automatic array.
	char *lastName;		// dynamic array, must allocate storage.

	cout << "Please enter the following information:" << endl;

	cout << "first name: ";
	cin >> firstName;

	lastName = new char[100];  // 'new' command allocates memory.
	cout << "last name: ";
	cin >> lastName;

	cout << "age: ";
	cin >> age;

	cout << "shoe size: ";
	cin >> shoeSize;

	bool canVote;  // can declare variables anywhere inside block.
	
	if (age - shoeSize >= 11)
		canVote = true;
	else 
		canVote = false;

	cout << firstName << ' ' << lastName;
	cout << ": according to our calculations you ";

	// hey, clever use of '?', who'd have thunk it !!.
	canVote ? cout << "can" : cout << "can't";

	cout << " vote. " << endl;

	// return dynamic storage to the system heap.
	delete [] lastName;

	return 0;
}
