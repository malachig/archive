/***
*	Main Program for Assignment 1.
*	Files: student_main.cpp, student_meth.cpp, string.cpp, string1.h, student.h
*	Author: Malachi Griffith
*	Date: Feb 2, 2002.
***/


#include <iostream>
using namespace std;

#include "student.h"  // Unless using namespaces, use 'traditional' #include
                      // <file_name.h> looks for in compiler set directories
                      // "file_name.h" looks for in *current* directory.

// Global Constants.
const int MAX_NUM_STUDENTS = 100; // At most 100 students.
const int MAX_STR_LENGTH = 100; // Max length of strings utilized.

// Menu choices available to user.
enum MenuChoiceEnum {ADD_STUDENT=1, REMOVE_STUDENT=2, PRINT_STUDENT=3, 
					 PRINT_LIST=4, EXIT=5};


// *** FUNCTION PROTOTYPES *** //

// Prints menu options, returns user selection.
MenuChoiceEnum get_menu_choice();

// Carries out menu action entered by user.
void perform_menu_action(MenuChoiceEnum menuChoice, Student *list[],
						 int &size);

// Functions for dealing with student objects.
Student* create_student();
int find_student(Student* list[], int &size);
bool add_student(Student *s, Student* list[], int &size, int maxSize);
void remove_student(Student* list[], int &size);
void print_student(Student* list[], int &size);
void print_list(Student* list[], int &size);

int 
main()
{
	Student *mList[MAX_NUM_STUDENTS]; // Array of pointers to students.
	MenuChoiceEnum  menuChoice;		  // User defined enumerated type var.
	
	// Size of array, will be passed to all other functions as a reference.
	int numStudents = 0;			  
									 
	// Continue to ask user for a selection until they select 'EXIT'.
	do
	{
		menuChoice = get_menu_choice();
		perform_menu_action( menuChoice, mList, numStudents );
	}
	while (menuChoice != EXIT);
  
	// Now delete any remaining students in the list.
	for (int i = 0; i < numStudents; i++)
		delete mList[i];

	numStudents = 0;  // Just a good practice.
	
	return 0;
}


/*** 
*	Prints menu options, returns user input.
*	Checks each entry for validity, repeats loop if the entry is
*	invalid.  This function accepts no arguments and returns a
*	MenuChoiceEnum value to main, where it will be sent to the function,
*	perform_menu_action().
***/

MenuChoiceEnum 
get_menu_choice()
{
  MenuChoiceEnum menuEntryNum = EXIT;  // Initialize choice to 'EXIT'
  char  ch; // If use int type for input, cin chokes on a letter.	

  cout << "\n\n";	
  do
  {
    cout << "Enter choice (1 - 5 ): \n";
    cout << "1 - ADD STUDENT\n"; 
    cout << "2 - REMOVE STUDENT\n"; 
    cout << "3 - PRINT STUDENT\n";
	cout << "4 - PRINT LIST\n";
    cout << "5 - EXIT\n" << "> "; 

    cin >> ch;			// Can't put into user defined type (yet)
    ch = ch - '0';	    // Convert ascii 'number' into digital number

    menuEntryNum = (MenuChoiceEnum)ch; // Assign to the user defined type.

  } // Repeat until valid entry.
  while ( menuEntryNum<ADD_STUDENT || menuEntryNum>EXIT );	

  return menuEntryNum;
}


/***
*	As the name suggests, carries out a menu action.  The menu entry number
*	from the user is accepted by this function and it carries out the action
*	requested by the user.
***/

void 
perform_menu_action(MenuChoiceEnum menuChoice, Student *list[], int &size)
{
	Student temp_student;

  switch( menuChoice )
  {
	case ADD_STUDENT:
		bool test;
	
		Student *tempStudent;
		tempStudent = create_student();	
    
		test = add_student(tempStudent, list, size, MAX_NUM_STUDENTS); 
		
		if (test == false)
		{
			// Student not added.
			cerr << "\n*** Student Not Added to List ***" << endl;
			cerr << "*** Array is Full or Student Pointer is Null ***\n";
		}

		break;
	
	case REMOVE_STUDENT:
		remove_student(list, size);
		break;

	case PRINT_STUDENT:
		print_student(list, size);
		break;

	case PRINT_LIST:
		print_list(list, size);
		break;
  }
}


// ** THE FOLLOWING FUNCTIONS MANIPULATE THE MASTER STUDENT LIST ** //

/***
*	The function create_student asks the user for student information,
*	It creates a Student object off the heap (using 'new') and sets
*	the appropriate data members and returns a pointer to the new student.
*	This student will then be added to the list by the add_student() function.
***/

Student* 
create_student()
{
	// Create a student object using Dynamic Memory Allocation.
	Student *newStudent;		// Student pointer.
	newStudent = new Student(); // Dynamically allocate memory for object.

	// If this allocation fails, new() returns null, skip rest of block:
	if (newStudent == 0)
		return newStudent;

	// Declare some automatic variables to deal with the input data. 
	char tempName[MAX_STR_LENGTH];	// Temporary student name storage.
	char tempAddr[MAX_STR_LENGTH];	// Temporary street address storage.
	char junk[1];					// Helps deal with input buffer issues.
	unsigned long tempNumber;		// Temporay student # storage.

	// Get the student data from the user.
	cout << "\nEnter the student name > ";
	cin.getline(junk, 1);  // Remove character already existing in buffer.
	cin.getline(tempName, MAX_STR_LENGTH);

	cout << "Enter the student's address > ";
	cin.getline(tempAddr, MAX_STR_LENGTH);
	
	cout << "Enter the student number > ";
	cin >> tempNumber;

	// Populate the Student Object using mutator functions.
	newStudent->set_student_name(tempName);
	newStudent->set_student_address(tempAddr);
	newStudent->set_student_number(tempNumber);

	return newStudent;
}


/***
*	If the Student pointer is not NULL, and the array of students is
*	not full, add the student to the end of the pointer array,
*	increment the 'size' argument accordingly.  Return 'true' is the
*   student was added, otherwise return false and print an error message.
*	Pre: A student object s, the student list[], the size of the array,
*		 and the max size allowed are all defined.  Size will be updated
*		 by reference, and the master list is also updated.
*	Post: A boolean value (true or false is returned).
***/

bool 
add_student(Student *s, Student* list[], int &size, int maxSize)
{
	if (s != NULL && size <= maxSize)
	{
		list[size] = s;  // Add student to the master list.
		size++;			 // Adjust the size of the master list.
		return true;	 // Student added successfully.
	}
	else 
		return false;	// Student not added.
}


/***
*	Asks the user for a student number, searches the array for a 
*	student with that number and returns the index of the array that
*	matches that student number.  Returns -1 is no student of that 
*	number is found.
*	Pre: The student list[] (array of pointers) and the size reference
*		 are defined.
*	Post: An integer value is returned.
***/

int 
find_student(Student* list[], int &size)
{
	int queryIndex = -1; // Location of target student in the master list.
	unsigned long studentNumQuery; // Student # to search for in list.
	unsigned long tempStudentNumber; // Temp student # for comparison.
	int i;  // loop control variable.

	cout << "\nEnter the student # > ";
	cin >> studentNumQuery;
	
	for (i = 0; i < size; i++)
	{
		// Get a student # from list and compare to query student #.
		tempStudentNumber = list[i]->get_student_number();

		if (studentNumQuery == tempStudentNumber)
		{
			queryIndex = i;
			break;  // exit loop if student is found.
		}
	}

	return queryIndex;  // Index location of target student.
}


/***
*	First calls the find_student() function to find the student object
*   to be targeted for removal.  If the student is found that student object 
*	is deleted from the array.  This leaves an empty space in the array so
*   a for loop is used to shift all of the students following the space so
*   that the array no longer contains any empty elements internally.  If the 
*   student was not found an error message to that effect is displayed.
*	Pre: The student list[] (array of pointers) and the size reference
*		 are defined.
***/

void 
remove_student(Student* list[], int &size)
{
	int indexValue;  // Index location of student to be removed.
	int i;			 // Loop control variable.

	indexValue = find_student(list, size);

	if (indexValue != -1)  // ie. if the student was found!
	{
		delete list[indexValue]; // Delete memory for that object.

		// Fill in the empty space created by shifting values "up".
		for (i = indexValue; i <= size; i++)
			list[i] = list[i+1];
		size--;
	}
	else
		cerr << "\n*** No Student of that Number in the List ***" << endl;
}


/***
*	First calls the find_student() function to find the student object
*   to be displayed to the user.  If it is found, the student is printed
*   using the print() behaviour of the Student Class.  If the student #
*	is not found, an error message to that effect is displayed.
*	Pre: The student list[] (array of pointers) and the size reference
*		 are defined.
***/

void 
print_student(Student* list[], int &size)
{
	int indexValue; // Array location of student to be printed.
	
	indexValue = find_student(list, size);

	// Now print the student using the print behaviour of the student object.
	if (indexValue != -1)  // ie. if the student was found!
		list[indexValue]->print();
	else 
		cerr << "\n*** No Student of that Number in the List ***" << endl;
}


/***
*	This function simply prints out the entire array of student objects from 
*	beginnning to end.  This ability is very helpful for testing the program. 
*   It provides a means of validating that a student has been successfully
*	removed.  It will also be helpful when sorting or other more involved
*	functions are added to the program (ie. for future feature developement).
*	If the list is empty, the user is informed.
*	Pre: The student list[] (array of pointers) and the size reference
*		 are defined.
***/

void 
print_list(Student* list[], int &size)
{
	int i;  // Loop control variable.

	for (i = 0; i < size; i++)
		list[i]->print();  // Print behaviour of student object.

	if (size == 0)
		cerr << "\n*** The list is Empty ***" << endl;
}
