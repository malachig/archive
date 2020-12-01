/***
*	Main Program for Assignment 2.
*	Files: main.cpp, student.cpp, course.cpp, my_string.cpp, linked_list.cpp
*	       student.h, course.h, my_string.h, linked_list.h
*	Author: Malachi Griffith
*	Date: March 11, 2002.
***/

#include <iostream>
#include <fstream>
#include <list>
using namespace std;

/***  I have included all of the classes being used in this file.  Even though
*     the course and linked list classes will be included within the student class
*     it is safest to include them all here in case this situation changes.  The 
*     class itself is set up to prevent multiple inclusions so that you don't have 
*     to worry about it.
***/
#include "student.h"  
#include "course.h"

// Global Constants and Variables
static list<Student*> g_StudentList;  // Creates the default empty student list
static list<Course*> g_CourseList;	  // Creates the default empty course list.
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
list<Student*>& master_student_list();
Student* create_student();
list<Student*>::iterator find_student(list<Student*> &sList);
void add_student();
void remove_student();
void print_student();
void print_list();

// Functions for dealing with courses.
list<Course*>& master_course_list();
void init_courses(list<Course*> &cList, String &fileName);
list<Course*>::iterator find_course(list<Course*> &cList);
void add_course(list<Student*> &sList, list<Course*> &cList);
void remove_course(list<Student*> &sList);
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
	list<Course*> &cList = master_course_list();  // Get reference to empty list.
	init_courses(cList, fileN);

	// Continue to ask the user for a selection until they select 'EXIT'.
	do
	{
		menuChoice = get_menu_choice();
		perform_menu_action(menuChoice);
	}
	while (menuChoice != EXIT);
  
	
	// Dynamic MEMORY CLEANUP for Student List. //
	
	list<Student*> &sList = master_student_list();  // Local reference to student list.

	if (!sList.empty())  // Only if the list is not empty.  Otherwise crash may occur. 
	{
		// Display number of students still in the list.
		cout << endl << sList.size() << " Students Still in Student List" << endl;
		cout << "Cleaning up Dynamic Memory" << endl;
		
		list<Student*>::iterator itr = sList.begin();

		while (itr != sList.end())
		{
			studentPtr = *itr;
			delete studentPtr;
			itr++;
		}
		// Now that objects are actually deleted, clear the list.
		sList.clear();
	}
	
	// Dynamic MEMORY CLEANUP for Course List. //
	cList = master_course_list();  // Reference to course list from above.

	if (!cList.empty())  // If not already empty.
	{
		cout << endl << cList.size() << " Courses Still in Course List" << endl;
		cout << "Cleaning up Dynamic Memory" << endl;

		list<Course*>::iterator itr = cList.begin();

		while (itr != cList.end())
		{
			coursePtr = *itr;
			delete coursePtr;
			itr++;
		}
		cList.clear();
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
	list<Student*> &sList = master_student_list();
	list<Course*> &cList = master_course_list();
	
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

list<Student*>& master_student_list()
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
	list<Student*> &sList = master_student_list();

	// Add student to the master list.
	sList.push_back(tempStudent);
}


/***
*	find_student():
*	Asks the user for a student number, searches the array for a 
*	student with that number and returns the index of the array that
*	matches that student number.  Returns null(0) if no student of that 
*	number is found.
*	Pre: A reference to the master list is passed is defined.
*	Post: An interator value is returned.  This specifies the location (by
*		  memory address) of the target student found.
***/

list<Student*>::iterator 
find_student(list<Student*> &sList)
{
	unsigned long studentNumQuery; // Student # to search for in list.
	unsigned long tempStudentNumber; // Temp student # for comparison.
	list<Student*>::iterator itr = sList.begin();

	Student *studentPtr;  // local pointer to student object

	cout << "\nEnter the student # > ";
	cin >> studentNumQuery;
	
	while (itr != sList.end())
	{
		studentPtr = *itr; 

		// Get a student # from list and compare to query student #.
		tempStudentNumber = studentPtr->get_student_number();

		if (studentNumQuery == tempStudentNumber)
			return itr;
		
		itr++;  // If not found, advance to next element in the list.
	}

	return 0;  // Return null if student not found.
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
	list<Student*>::iterator studentFound; // Location of student to be removed.
	list<Student*>::iterator eraseStudent; // For use with list.erase member function.
	bool test;		 // Test variable for removal of students from course lists.
	
	Student *tempStudent;  // Points to any student object
	Course *tempCourse;	   // Points to any course object
	list<Student*> &sList = master_student_list(); // reference to master s list.
	list<Course*> &cList = master_course_list();   // reference to master c list. 

	if (sList.empty())  // Attempts to acces an empty list may cause crashes.
	{
		cerr << "\n*** The List is Empty ***" << endl;
		return;
	}

	// Ask the user for a student number and search for that student.
	// Get the iterator value (address) for that element in the list.
	studentFound = find_student(sList);

	if (studentFound == 0)  // ie. if the student was not found!
	{
		cerr << "\n*** No Student of that Number in the List ***" << endl;
		return;
	}

	// If the student WAS found then do the following.
	tempStudent = *studentFound;  // Get pointer to that student using iterator

	
	// Remove this student from all course lists.

	list<Course*>::iterator itr = cList.begin();

	while (itr != cList.end())
	{
		// Operator overloading, '*' actually returns an element in list
		tempCourse = *itr;

		// Remove the Student from the Course's List of Students
		// Note that remove_student checks if the student is in a particular
		// course before doing anything.
		test = tempCourse->remove_student(tempStudent);
		
		itr++;  
	}

	// Remove student from the list and get iterator of next student in list.
	eraseStudent = sList.erase(studentFound);

	// Now that the course lists have been cleaned up we can delete the student.
	// *** Hopefully the erase method does not delete memory itself ***
	delete tempStudent;
		
	//  NO NEED to fill in the empty space created by shifting values "up".
	//  The STL handles all this for us!
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
	Student *studentPtr;  // local pointer to student object
	list<Student*>::iterator studentFound; // Location of student to be printed
	
	list<Student*> &sList = master_student_list();  //local reference to master list.

	if (sList.empty())  // Attempts to acces an empty list may cause crashes.
	{
		cerr << "\n*** The List is Empty ***" << endl;
		return;
	}

	studentFound = find_student(sList);

	// Now print the student using the print behaviour of the student object.
	if (studentFound == 0)  // ie. if the student was not found!
	{
		cerr << "\n*** No Student of that Number in the List ***" << endl;
		return;
	}
	
	// Get the student object from the list
	studentPtr = *studentFound;  // Get pointer to that student using iterator

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
	list<Student*> &sList = master_student_list();  //local reference to student list.
	list<Course*> &cList = master_course_list();	  //local reference to course list.
	Student *studentPtr;  // local pointer to student object
	Course *coursePtr;	  // local pointer to course object
	
	if (sList.empty())  // Attempts to acces an empty list are avoided.
		cerr << "\n*** No Students in the List ***\n" << endl;
	
	list<Student*>::iterator itr = sList.begin();

	if (!sList.empty())
	{
		cout << "\nSTUDENTS CURRENTLY IN THE DATABASE" << endl;
		// Print out the student objects.
		while (itr != sList.end())
		{
			// Get a student object from the list
			studentPtr = *itr; // uses operator overloading for '[]'
	
			// Now print the student using the print behaviour of the student object.
			studentPtr->print();
			itr++;
		}
	}

	if (cList.empty())  // Attempts to acces an empty list are avoided.
		cerr << "\n*** No Courses in the List ***" << endl;

	list<Course*>::iterator itr_c = cList.begin();

	if (!cList.empty())
	{
		cout << "ENROLLMENT FOR ALL COURSES IN THE DATABASE\n" << endl;
		// Print out the course objects.
		while (itr_c != cList.end())
		{
			// Get a course object from the list
			coursePtr = *itr_c;

			// Now print the course using the print behaviour of the course object.
			coursePtr->print();
			itr_c++;
		}
	}
}


// ** THE FOLLOWING FUNCTIONS MANIPULATE THE MASTER COURSE LIST ** //

/***
*	This function returns a reference to the global master course list.
*	Same idea as an earlier function for student list.
***/
list<Course*>& master_course_list()
{
	return g_CourseList;
}

/***
*	init_courses():
*	Accepts a reference to a course list and populates it by accessing the 
*	file specified by the user.  Initializes the list of courses and adds them 
*	to the given course list.
*	Pre: A reference to the master course list is defined as well as a filename string.
***/
void 
init_courses(list<Course*> &cList, String &fileName)
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
		// Add the course to the STL list
		cList.push_back(c);
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

list<Course*>::iterator
find_course(list<Course*> &cList)
{
	int counter = 1;
	int index = -1;	// Index of course selected. (initialize to -1)
	Course *c;		// Local pointer to course object.


	// First display the courses available to the user:
	cout << "\n\nCourses Available: \n\n";

	// Use iterator to 'loop over the list.
	list<Course*>::iterator itr = cList.begin();

	while (itr != cList.end())
	{
		c = *itr;
		cout << counter << ". - " << c->get_name() << endl;
		counter++;
		itr++;
	}

	cout << "Select a Course > ";
	cin >> index;

	// Make sure their entry was valid, otherwise return -1
	if (index <= 0 || index > cList.size())  // ie. not 1-5.
	{
		cerr << "\n*** Not a Valid Course Entry ***" << endl;
		return 0;
	}

	// Now loop through the list again and grab the iterator corresponding to
	// the user's selection.
	itr = cList.begin();
	counter = 1;

	while (itr != cList.end())
	{
		if (counter == index)
			return itr;
		counter++;
		itr++;
	}

	return 0; // Return Null if course not found.
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
add_course(list<Student*> &sList, list<Course*> &cList)
{
	list<Course*>::iterator courseFound; // Location of course.
	list<Student*>::iterator studentFound;  // Location of student.
	Student *tempStudent; // Pointer to target student object.
	Course *tempCourse;	  // Pointer to target course object.
	bool test;			  // Bool variable, for testing whether student was 
						  // successfully added to a Course's List.

	// If the student list is empty, no point in continuing.
	if (sList.empty())
	{
		cerr << "\n*** No Students in the List ***" << endl;
		return;
	}

	// Make sure there are courses in the list.
	if (cList.empty())
	{
		cerr << "\n*** No Courses in the List ***" << endl;
		return;
	}

	// First get the student selection.
	studentFound = find_student(sList);

	// If that student was not found, display error message and return.
	if (studentFound == 0)
	{
		cerr << "\n*** No Student of that Number in the List ***" << endl;
		return;
	}

	tempStudent = *studentFound;
	
	// Get the course selection.
	courseFound = find_course(cList);

	// If that course was not found, return.
	if (courseFound == 0)
		return;
	
	tempCourse = *courseFound;

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
remove_course(list<Student*> &sList)
{
	list<Course*>::iterator courseFound; // Location of course.
	list<Student*>::iterator studentFound;  // Location of student.
	Student *tempStudent; // Pointer to target student object.
	Course *tempCourse;	  // Pointer to target course object.
	bool test;			  // Bool variable, for testing whether student was 
						  // successfully removed from a Course's List.
	
	// Get the student of interest from the user.
	studentFound = find_student(sList);
	
	// If the student was not found display error message and return.
	if (studentFound == 0)
	{
		cerr << "\n*** No Student of that Number in the List ***" << endl;
		return;
	}

	tempStudent = *studentFound;

	// Get a local copy of the student's course list
	list<Course*> &cList = tempStudent->get_course_list();

	// Check if course list is empty for this student, if so, return.
	if (cList.empty())
	{
		cerr << "\n*** Student is not enrolled in any courses ***" << endl;
		return;
	}

	// Get the course selection.
	courseFound = find_course(cList);

	// If that course was not found, return.
	if (courseFound == 0)
		return;
	
	tempCourse = *courseFound;


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
	list<Course*>::iterator courseFound; // Location of course.
	Course *tempCourse;	  // Pointer to target course object.
	
	list<Course*> &cList = master_course_list();	  //local reference to course list.
	
	if(cList.empty())
	{
		cerr << "\n*** The Course List is Empty ***" << endl;
		return;
	}

	// Ask user which Course they want to see the enrollment for.
	courseFound = find_course(cList);

	// Check for error in find_course result.
	if (courseFound == 0)
		return;  // Error message already displayed in find_course()

	tempCourse = *courseFound;

	// Display the Course's student list using the method designed to do just that.
	tempCourse->print_students();
}






