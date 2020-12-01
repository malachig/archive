/***
*	Main Program for Assignment 2.
*	Files: main.cpp, student.cpp, course.cpp, my_string.cpp, linked_list.cpp
*	       student.h, course.h, my_string.h, linked_list.h
*	Author: Malachi Griffith
*	Date: March 11, 2002.
***/

#include <iostream>
using namespace std;
#include <fstream>
using namespace std;

/***  I have included all of the classes being used in this file.  Even though
*     the course and linked list classes will be included within the student class
*     it is safest to include them all here in case this situation changes.  The 
*     class itself is set up to prevent multiple inclusions so that you don't have 
*     to worry about it.
***/
#include "linked_list.h"
#include "student.h"  
#include "course.h"

// Global Constants and Variables
static LinkedList g_StudentList;  // Creates the default empty student list
static LinkedList g_CourseList;	  // Creates the default empty course list.
const int MAX_STR_LENGTH = 100; // Max length of strings utilized.

// Menu choices available to user.
enum MenuChoiceEnum {ADD_STUDENT=1, REMOVE_STUDENT=2, ADD_COURSE=3, 
					 REMOVE_COURSE=4, PRINT_STUDENT=5, PRINT_COURSE=6, 
					 PRINT_LIST=7, EXIT=8};


// *** FUNCTION PROTOTYPES *** //

// Prints menu options, returns user selection.
MenuChoiceEnum get_menu_choice();

// Carries out menu action entered by user.
void perform_menu_action(MenuChoiceEnum menuChoice);

// Functions for dealing with students.
LinkedList & master_student_list();
Student* create_student();
int find_student(const LinkedList &sList);
void add_student();
void remove_student();
void print_student();
void print_list();

// Functions for dealing with courses.
LinkedList & master_course_list();
void init_courses(LinkedList &cList, String &fileName);
int find_course(const LinkedList &cList);
void add_course(LinkedList &sList, LinkedList &cList);
void remove_course(LinkedList &sList);
void print_course();

// *** MAIN *** //

int 
main(int argc, char *argv[])
{
	MenuChoiceEnum  menuChoice;	// User defined enumerated type variable
	Student *studentPtr;		// Points to any student object  
	Course *coursePtr;			// Points to any course object
	char tempFileName[100];
	String fileN;

	fileN = argv[1];
	
	// Check for lack of input arguments (ie. the neccessary filename)
	if (argc < 2)
	{
		cerr << "\n*** Need an Input File Containing Course Info ***" << endl;
		cout << "Enter the filename > ";
		cin >> tempFileName;
		fileN = tempFileName;
	}

	// Initialize the course list.
	LinkedList &cList = master_course_list();  // Get reference to empty list.
	init_courses(cList, fileN);

	// Continue to ask the user for a selection until they select 'EXIT'.
	do
	{
		menuChoice = get_menu_choice();
		perform_menu_action(menuChoice);
	}
	while (menuChoice != EXIT);
  
	
	// Dynamic MEMORY CLEANUP for Student List. //
	
	LinkedList &sList = master_student_list();  // Local reference to student list.

	if (sList.size() > 0)  // Only if the list is not empty.  Otherwise crash may occur. 
	{
		// Display number of students still in the list.
		cout << endl << sList.size() << " Students Still in Student List" << endl;
		cout << "Cleaning up Dynamic Memory" << endl;
		
		for (int i = 0; i < sList.size(); i++)
		{
			studentPtr = (Student*)sList[i]; // uses operator overloading for '[]'
			delete studentPtr;
		}
		// Clear list to remove the nodes that used to point to the objects deleted.
		sList.clear_list();
	}
	
	// Dynamic MEMORY CLEANUP for Course List. //
	cList = master_course_list();  // Reference to course list from above.

	if (cList.size() > 0)  // If not already empty.
	{
		cout << endl << cList.size() << " Courses Still in Course List" << endl;
		cout << "Cleaning up Dynamic Memory" << endl;

		for (int i = 0; i < cList.size(); i++)
		{
			coursePtr = (Course*)cList[i];
			delete coursePtr;
		}
		// Clear list to remove nodes.
		cList.clear_list();
	}
	return 0;
}


// *** GIVING THE USER A MENU CHOICE AND PERFORMING THAT ACTION  *** //

/*** 
*	get_menu_choice():
*	Prints menu options, returns user input.
*	Checks each entry for validity, repeats loop if the entry is
*	invalid.  This function accepts no arguments and returns a
*	MenuChoiceEnum value to main, where it will be sent to the function,
*	perform_menu_action().
*	Post: The enumerated type MenuChoice value is returned.
***/

MenuChoiceEnum 
get_menu_choice()
{
  MenuChoiceEnum menuEntryNum = EXIT;  // Initialize choice to 'EXIT'
  char  ch; // If use int type for input, cin chokes on a letter.	

  cout << "\n\n\n\n";	
  do
  {
    cout << "Enter choice (1 - 5): \n";
    cout << "1 - ADD STUDENT\n"; 
    cout << "2 - REMOVE STUDENT\n"; 
	cout << "3 - ADD COURSE\n";
	cout << "4 - REMOVE COURSE\n";
    cout << "5 - PRINT STUDENT\n";
	cout << "6 - PRINT COURSE ENROLLMENT\n";
	cout << "7 - PRINT LIST\n";
    cout << "8 - EXIT\n" << "> "; 

    cin >> ch;			// Can't put into user defined type (yet)
    ch = ch - '0';	    // Convert ascii 'number' into digital number

    menuEntryNum = (MenuChoiceEnum)ch; // Assign to the user defined type.

  } // Repeat until valid entry.
  while (menuEntryNum<ADD_STUDENT || menuEntryNum>EXIT);	

  return menuEntryNum;
}


/***
*	perform_menu_action():
*	As the name suggests, carries out a menu action.  The menu entry number
*	from the user is accepted by this function and it carries out the action
*	requested by the user.
*	Pre: The menuChoice value chosen by the user is defined.
***/

void 
perform_menu_action(MenuChoiceEnum menuChoice)
{
	LinkedList &sList = master_student_list();
	LinkedList &cList = master_course_list();
	
  switch( menuChoice )
  {
	case ADD_STUDENT:
		// Simply call add_student() function which will in turn call 
		// create_student() and add the new student to the linked list.
		add_student();
		break;
		
	case REMOVE_STUDENT:
		remove_student();
		break;
	
	case ADD_COURSE:
		add_course(sList, cList);
		break;
	
	case REMOVE_COURSE:
		remove_course(sList);
		break;

	case PRINT_COURSE:
		print_course();
		break;
	
	case PRINT_STUDENT:
		print_student();
		break;
	
	case PRINT_LIST:
		print_list();
		break;
  }
}


// ** THE FOLLOWING FUNCTIONS MANIPULATE THE MASTER STUDENT LIST ** //


/***
*	master_student_list():
*	This function returns a reference to the global master student list.
*	Can be used to create a reference to the list which can passed to  
*	functions that need to access the list.  Or it can be used in the 
*	function call itself or even within the function being called.
*	Post: Returns a reference to a LinkedList object.
***/

LinkedList & master_student_list()
{
	return g_StudentList;
}


/***
*	create_student():
*	The function create_student() asks the user for student information,
*	It creates a Student object off the heap (using 'new') and sets
*	the appropriate data members and then returns a pointer to the new student.
*	This student will then be added to the list by the add_student() function.
***/

Student* 
create_student()
{
	// Declare some automatic variables to deal with the input data. 
	char tempName[MAX_STR_LENGTH];	// Temporary student name storage.
	char tempAddr[MAX_STR_LENGTH];	// Temporary street address storage.
	char tempAdvisor[MAX_STR_LENGTH];  // Temp advisor name storage.
	char tempThesis[MAX_STR_LENGTH];   // Temp thesis name storage.
	char tempProject[MAX_STR_LENGTH];  // Temp project name storage.
	char junk[1];					// Helps deal with input buffer issues.
	unsigned long tempNumber;		// Temporay student # storage.
	int student_type = 0;		    // 1 = GradStudent, 2 = UnderGradStudent.

	Student *newStudent;		// Student pointer.

	// Get the student data from the user.
	cout << "\nEnter the student name > ";
	cin.getline(junk, 1);  // Remove character already existing in buffer.
	cin.getline(tempName, MAX_STR_LENGTH-1);

	cout << "Enter the student's address > ";
	cin.getline(tempAddr, MAX_STR_LENGTH-1);
	
	cout << "Enter the student number > ";
	cin >> tempNumber;

	while (student_type > 2 || student_type < 1)
	{
		cout << "Is the student a GradStudent(1) or UnderGrad Student (2)? > ";
		cin >> student_type;
		if (student_type > 2 || student_type < 1)
			cerr << "\n*** Invalid Selection ***" << endl;
	}

	if (student_type == 1)
	{
		cout << "Enter the Grad's advisor name > ";
		cin.getline(junk, 1);  // Remove character already existing in buffer.
		cin.getline(tempAdvisor, MAX_STR_LENGTH-1);
		cout << "Enter the Grad's thesis title > ";
		cin.getline(tempThesis, MAX_STR_LENGTH-1);
			
		// Create a GRAD student object using Dynamic Memory Allocation.
		newStudent = new GradStudent(tempThesis, tempAdvisor);

		// If this allocation fails, new() returns null, skip rest of block:
		if (newStudent == 0)
		{
			cerr << "\n*** Memory allocation for student failed ***" << endl;
			return newStudent;
		}

		// Assign Person Attributes gathered above
		newStudent->set_name(tempName);
		newStudent->set_address(tempAddr);
		newStudent->set_student_number(tempNumber);
	}

	if (student_type == 2)
	{
		cout << "Enter the UnderGrad's project name > ";
		cin.getline(junk, 1);  // Remove character already existing in buffer.
		cin.getline(tempProject, MAX_STR_LENGTH-1);
				
		// Create an UNDERGRAD student object using Dynamic Memory Allocation.
		newStudent = new UnderGradStudent(tempProject);

		// If this allocation fails, new() returns null, skip rest of block:
		if (newStudent == 0)
		{
			cerr << "\n*** Memory allocation for student failed ***" << endl;
			return newStudent;
		}

		// Assign Person Attributes gathered above
		newStudent->set_name(tempName);
		newStudent->set_address(tempAddr);
		newStudent->set_student_number(tempNumber);
	}

	return newStudent;
}


/***
*	add_student():
*	If the Student pointer is not NULL add the student to the master linked 
*	list of students by first getting a reference to that list.  
*	Students can be added to this list indefinitely.  Since it is a linked list,
*	in theory it could grow until it used up all available memory.
***/

void
add_student()
{
	Student *tempStudent;  // Pointer to student object.

	// Create the student.
	tempStudent = create_student();  

	// Ensure that the student was created successfully.
	if (tempStudent == NULL)
		return;
	
	// Get local access to the master list via a reference
	LinkedList &sList = master_student_list();

	// Add student to the master list.
	sList.add(tempStudent);
}


/***
*	find_student():
*	Asks the user for a student number, searches the array for a 
*	student with that number and returns the index of the array that
*	matches that student number.  Returns -1 is no student of that 
*	number is found.
*	Pre: A reference to the master list is passed is defined.  Since this
*	     function is only searching the list and not changing values, the 
*		 reference is made 'const' to prevent accidental changes to the list.
*	Post: An integer value is returned.  This integer specifies the index of the
*		  target student found.
***/

int 
find_student(const LinkedList &sList)
{
	int queryIndex = -1; // Location of target student in the master list.
	unsigned long studentNumQuery; // Student # to search for in list.
	unsigned long tempStudentNumber; // Temp student # for comparison.
	int i;  // loop control variable.
	
	Student *studentPtr;  // local pointer to student object

	cout << "\nEnter the student # > ";
	cin >> studentNumQuery;
	
	for (i = 0; i < sList.size(); i++)
	{
		studentPtr = (Student*)sList[i]; // uses operator overloading for '[]'

		// Get a student # from list and compare to query student #.
		tempStudentNumber = studentPtr->get_student_number();

		if (studentNumQuery == tempStudentNumber)
		{
			queryIndex = i;
			break;  // exit loop if student is found.
		}
	}

	return queryIndex;  // Index location of target student.
}


/***
*	remove_student():
*	First calls the find_student() function to find the student object
*   to be targeted for removal.  If the student is found, that student object 
*	is deleted from the array.  Since the objects are tracked as a linked list,
*   no organization of the remaining students is required after removing one.  
*   If the student was not found an error message to that effect is displayed.
***/

void 
remove_student()
{
	int indexValue;  // Index location of student to be removed.
	bool test;		 // Test variable for removal of students from course lists.
	
	Student *tempStudent;  // Points to any student object
	Course *tempCourse;	   // Points to any course object
	LinkedList &sList = master_student_list();  // reference to master student list.
	LinkedList &cList = master_course_list();   // reference to master course list. 

	if (sList.size() <= 0)  // Attempts to acces an empty list may cause crashes.
	{
		cerr << "\n*** The List is Empty ***" << endl;
		return;
	}

	// Ask the user for a student number and search for that student.
	indexValue = find_student(sList);

	if (indexValue == -1)  // ie. if the student was not found!
	{
		cerr << "\n*** No Student of that Number in the List ***" << endl;
		return;
	}

	// If the student WAS found then do the following.
	
	// Remove student from the list and get a pointer to the student object.
	tempStudent = (Student*)sList.remove(indexValue);
	
	// Remove this student from all course lists.
	for (int i = 0; i < cList.size(); i++)
	{
		tempCourse = (Course*)cList[i];

		// Remove the Student from the Course's List of Students
		// Note that remove_student checks if the student is in a particular
		// course before doing anything.
		test = tempCourse->remove_student(tempStudent);
	}

	// Now that the course lists have been cleaned up we can delete the student.
	delete tempStudent;
		
	//  NO NEED to fill in the empty space created by shifting values "up".
	//  Since this is a linked list the objects are not placed sequentially
	//  in memory.  The remove function just rearranges the links so that the 
	//  object removed is no longer part of the link.
}


/***
*	print_student():
*	First calls the find_student() function to find the student object
*   to be displayed to the user.  If it is found, the student is printed
*   using the print() behaviour of the Student Class.  If the student #
*	is not found, an error message to that effect is displayed.
***/

void 
print_student()
{
	int indexValue; // Index location of student to be printed.
	Student *studentPtr;  // local pointer to student object
	
	const LinkedList &sList = master_student_list();  //local reference to master list.

	if (sList.size() <= 0)  // Attempts to acces an empty list may cause crashes.
	{
		cerr << "\n*** The List is Empty ***" << endl;
		return;
	}

	indexValue = find_student(sList);

	// Now print the student using the print behaviour of the student object.
	if (indexValue == -1)  // ie. if the student was not found!
	{
		cerr << "\n*** No Student of that Number in the List ***" << endl;
		return;
	}
	
	// Get the student object from the list
	studentPtr = (Student*)sList[indexValue];

	// Print the student info using the student object print() behaviour.
	studentPtr->print();
}


/***
*	print_list():
*	This function simply prints out the entire list of student objects from 
*	beginnning to end.  This ability is very helpful for testing the program. 
*   It provides a means of validating that a student has been successfully
*	removed.  It will also be helpful when sorting or other more involved
*	functions are added to the program (ie. for future feature developement).
*	If the list is empty, the user is informed.
***/

void 
print_list()
{
	int i;
	const LinkedList &sList = master_student_list();  //local reference to student list.
	const LinkedList &cList = master_course_list();	  //local reference to course list.
	Student *studentPtr;  // local pointer to student object
	Course *coursePtr;	  // local pointer to course object
	
	if (sList.size() <= 0)  // Attempts to acces an empty list are avoided.
		cerr << "\n*** No Students in the List ***\n" << endl;
		
	if (sList.size() > 0)
	{
		cout << "\nSTUDENTS CURRENTLY IN THE DATABASE" << endl;
		// Print out the student objects.
		for (i = 0; i < sList.size(); i++)
		{
			// Get a student object from the list
			studentPtr = (Student*)sList[i]; // uses operator overloading for '[]'
	
			// Now print the student using the print behaviour of the student object.
			studentPtr->print();
		}
	}

	if (cList.size() <= 0)  // Attempts to acces an empty list are avoided.
		cerr << "\n*** No Courses in the List ***" << endl;

	if (cList.size() > 0)
	{
		cout << "ENROLLMENT FOR ALL COURSES IN THE DATABASE\n" << endl;
		// Print out the course objects.
		for (i = 0; i < cList.size(); i++)
		{
			// Get a course object from the list
			coursePtr = (Course*)cList[i];

			// Now print the course using the print behaviour of the course object.
			coursePtr->print();
		}
	}
}


// ** THE FOLLOWING FUNCTIONS MANIPULATE THE MASTER COURSE LIST ** //

/***
*	This function returns a reference to the global master course list.
*	Same idea as an earlier function for student list.
***/

LinkedList & master_course_list()
{
	return g_CourseList;
}

/***
*	init_courses():
*	Accepts a reference to a course list and populates it with courses by hard
*	coding for now.  Initializes the list of courses and adds them to the given
*	course list.
*	Pre: A reference to the master course list is defined.
***/
void 
init_courses(LinkedList &cList, String &fileName)
{
	int numEntries;
	char tempName[1000];
	char junk[5];
	const char *fileN;

	fileN = fileName.data();
	
	ifstream iFile(fileN);

	// Test if file was opened successfully
	if (!iFile.is_open())
	{
		cerr << "\n*** Course File Not Opened Properly ***" << endl;
		exit(0);
	}

	// Create 5 Course Pointers 
	Course *c;

	iFile >> numEntries;  // Get number from beginning of file;
	iFile.getline(junk, 4);  // Clear character left in buffer.

	for (int i = 0; i < numEntries; i++)
	{
		iFile.getline(tempName, 999);

		// Create Course objects via Dynamic Memory Allocation
		c = new Course(tempName);

		// If the memory allocation fails (new() returns null) then return.
		if (c == 0)
		{
			cerr << "\n*** Memory Allocation Failure ***" << endl;
			return;
		}
		// Add the course to the linked list
		cList.add(c);
	}
	// Close the input file
	iFile.close();
}


/***
*	find_course():
*	Displays the given list of courses available in the master course list,
*	asks the user to select one (by index # starting at 1) returns the 
*	selected *list* index, starting at 0.  Returns -1 if course not found,
*	ie. in this case if the user selects a number that is not in the list.
*	Pre: A reference to the master course list is defined.
*	Post: An integer value is returned.  This integer specifies the index of the
*		  target course found.
***/

int 
find_course(const LinkedList &cList)
{
	int index = -1;	// Index of course selected. (initialize to -1)
	Course *c;		// Local pointer to course object.


	// First display the courses available to the user:
	cout << "\n\nCourses Available: \n\n";

	for (int i = 0; i < cList.size(); i++)
	{
		c = (Course*)cList[i];
		cout << i+1 << ". - " << (*c) << endl;
	}

	cout << "Select a Course > ";
	cin >> index;

	// Make sure their entry was valid, otherwise return -1
	if (index <= 0 || index > cList.size())  // ie. not 1-5.
	{
		cerr << "\n*** Not a Valid Course Entry ***" << endl;
		return -1;
	}

	return (index - 1); // Adjust value because the user selected from 1-5!
}



/***
*	add_course():
*	Finds the student (specified by the user) in the given student 
*	list (sList).  If the student is found, attempts to add a course 
*	(also specified by the user) from the given course list (cList)
*	to the found student.  Calls the function find_course() to get a 
*	course and the function find_student() to get a student.
*	Pre: References to the master course and student lists are defined.
***/

void 
add_course(LinkedList &sList, LinkedList &cList)
{
	int sIndex;			  // Index location of target student.
	int cIndex;			  // Index location of target course.
	Student *tempStudent; // Pointer to target student object.
	Course *tempCourse;	  // Pointer to target course object.
	bool test;			  // Bool variable, for testing whether student was 
						  // successfully added to a Course's List.

	// If the student list is empty, no point in continuing.
	if (sList.size() == 0)
	{
		cerr << "\n*** No Students in the List ***" << endl;
		return;
	}

	// Make sure there are courses in the list.
	if (cList.size() <= 0)
	{
		cerr << "\n*** No Courses in the List ***" << endl;
		return;
	}

	// First get the student selection.
	sIndex = find_student(sList);

	// If that student was not found, display error message and return.
	if (sIndex == -1)
	{
		cerr << "\n*** No Student of that Number in the List ***" << endl;
		return;
	}

	tempStudent = (Student*)sList[sIndex];
	
	// Get the course selection.
	cIndex = find_course(cList);

	// If that course was not found, return.
	if (cIndex == -1)
		return;
	
	tempCourse = (Course*)cList[cIndex];

	// Since everything is okay, add the course to the student's personal list
	tempStudent->add_course(tempCourse);

	// Now add the student to the Course's List of Students
	test = tempCourse->add_student(tempStudent);

	// Error checking.
	if (test)
		cout << "\nEnrollment for this Course updated successfully" << endl;
	else
		cerr << "\n*** Update of Enrollment for this course failed ***" << endl;
}


/***
*	remove_course():
*	Finds the student in the given list, if student found, display list
*	of courses this student is enrolled in.  Then ask the user which course
*	to remove from the student's list and remove it.  Beware the gauntlet of 
*	error checking.
*	Pre: A reference to the master student list is defined.
***/

void 
remove_course(LinkedList &sList)
{
	int sIndex;			  // Index location of target student.
	int cIndex;			  // Index location of target course.
	Student *tempStudent; // Pointer to target student object.
	Course *tempCourse;	  // Pointer to target course object.
	bool test;			  // Bool variable, for testing whether student was 
						  // successfully removed from a Course's List.
	
	// Get the student of interest from the user.
	sIndex = find_student(sList);
	
	// If the student was not found display error message and return.
	if (sIndex == -1)
	{
		cerr << "\n*** No Student of that Number in the List ***" << endl;
		return;
	}

	tempStudent = (Student*)sList[sIndex];

	// Get a local copy of the student's course list
	LinkedList &cList = tempStudent->get_course_list();

	// Check for course list is empty for this student, if so, return.
	if (cList.size() == 0)
	{
		cerr << "\n*** Student is not enrolled in any courses ***" << endl;
		return;
	}

	// Display this student's courses.
	tempStudent->print_courses();

	// Ask user which course to remove
	cout << "Which course would you like to remove (enter number) > ";
	cin >> cIndex;

	// Check that the user's entry is valid
	if (cIndex > cList.size() || cIndex <= 0)
	{
		cerr << "\n*** Not a Valid Course Entry ***" << endl;
		return;
	}

	tempCourse = (Course*)cList[cIndex-1];  // Must adjust index value.

	// Now that everything has been checked, actually remove the course.
	// Do not delete memory for that course because student does not OWN 
	// the course objects.
	tempStudent->remove_course(tempCourse);

	// Remove the Student from the Course's List of Students
	test = tempCourse->remove_student(tempStudent);

	// Error checking.
	if (test)
		cout << "\nEnrollment for this Course updated successfully" << endl;
	else
		cerr << "\n*** Update of Enrollment for this course failed ***" << endl;
}

void print_course()
{
	int cIndex;	 // Index location of target course.
	Course *tempCourse;	  // Pointer to target course object.
	
	const LinkedList &cList = master_course_list();	  //local reference to course list.
	
	if(cList.size() < 1)
	{
		cerr << "\n*** The Course List is Empty ***" << endl;
		return;
	}

	// Ask user which Course they want to see the enrollment for.
	cIndex = find_course(cList);

	// Check for error in find_course result.
	if (cIndex == -1)
		return;  // Error message already displayed in find_course()

	tempCourse = (Course*)cList[cIndex];

	// Display the Course's student list using the method designed to do just that.
	tempCourse->print_students();
}






