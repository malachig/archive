#include <iostream>
#include <string>
#include <cstring>  // To use ci.getline
using namespace std;

#include "student.h" 
#include "course.h"


// global constants
const int MAX_NUM_STUDENTS = 10;  // at most 10 students for now.
const int MAX_NUM_COURSES = 10; // at most 10 courses for now.
const int MAX_LENGTH = 100; // max length of course name.

// global variables, try to minimize these.
// create and initialize an array of 10 student objects.
Student *studentListG = new Student[MAX_NUM_STUDENTS];  // Master student list.
int numStudentG = 0;  // Initial number of students in master list.

// create and initialize an array of 10 course objects.
Course *courseListG = new Course[MAX_NUM_COURSES];
int numCourseG = 0;

/*  Note to Self:  It would be better to create and initialize students
*   courses only as they are needed.  You don't know how many students 
*   you are going to have ultimately.
*/

// Menu choices available to user.
enum MenuChoiceEnum {ADD_STUDENT, ADD_COURSE, PRINT_STUDENT_LIST, 
		     PRINT_COURSE_LIST, QUIT};

//  *** Function prototypes ***

// Prints menu options, returns user input
MenuChoiceEnum get_menu_choice();

// Carries out menu action entered by user:
void perform_menu_action(MenuChoiceEnum menuChoice);
void add_student();
void add_course();
void print_student_list();
void print_course_list();


int main()
{
	MenuChoiceEnum menuChoice;

	do
	{
		menuChoice = get_menu_choice();
		perform_menu_action(menuChoice);
	}
	while (menuChoice != QUIT);

	return 0;
	
	// Free the dynamically assigned memory
	delete [] studentListG; 
	delete [] courseListG;
}

// Prints menu options, returns user input
MenuChoiceEnum get_menu_choice()
{
	MenuChoiceEnum menuEntryNum = QUIT;
	char ch;  // if use int type, cin chokes on a letter.

	cout << "\n\n";

	do
	{
		cout << "Enter choice (0 - 4):\n";
		cout << "0 - ADD STUDENT\n";
		cout << "1 - ADD COURSE\n";
		cout << "2 - PRINT STUDENT LIST\n";
		cout << "3 - PRINT COURSE LIST\n";
		cout << "4 - QUIT\n" << "> ";

		cin >> ch;  // can't put into user defined type (yet).
		cout << endl;

		ch = ch - '0';  // convert ascii 'number' into digital number.

		// Now assign it to the user defined type.
		menuEntryNum = (MenuChoiceEnum)ch;
	}
	// Repeat until the entry is valid.
	while (menuEntryNum < ADD_STUDENT || menuEntryNum > QUIT);

	return menuEntryNum;
}

// As the name suggests, carry out the menu action.
void perform_menu_action(MenuChoiceEnum menuChoice)
{
	switch(menuChoice)
	{
		case ADD_STUDENT:
			add_student();
			break;

		case ADD_COURSE:
			add_course();
			break;

		case PRINT_STUDENT_LIST:
			print_student_list();
			break;

		case PRINT_COURSE_LIST:
			print_course_list();
			break;
	}
}

// ** The following functions manipulate the master student list ** //

void add_student()
{
	char temp_name[MAX_LENGTH];
	char junk[2];
	int temp_number;

	cout << "Enter the name of the student > ";
	cin.getline(junk, 1);
	cin.getline(temp_name, MAX_LENGTH);

	cout << "Enter the student number > ";
	cin >> temp_number;
	
	studentListG[numStudentG].set_student_name(temp_name);
	studentListG[numStudentG].set_student_number(temp_number);	

	numStudentG++;

}

void add_course()
{
	char temp_course[MAX_LENGTH];
	char junk[2];

	double temp_number;

	cout << "Enter the full name of the course > ";
	cin.getline(junk, 1); // Advances the input cursor past endline
						  // character in input buffer already.

	cin.getline(temp_course, MAX_LENGTH);
	cout << "Enter the course number > ";
	cin >> temp_number;

	courseListG[numCourseG].set_course_name(temp_course);
	courseListG[numCourseG].set_course_number(temp_number);

	numCourseG++;
}

void print_student_list()
{
	int i;

	for (i = 0; i < numStudentG; i++)
	{
		cout << "Student " << i+1 << " is: " << studentListG[i].get_student_name() << endl;
		cout << "Student # is: " << studentListG[i].get_student_number() << endl << endl;
	}
}

void print_course_list()
{

	int i;

	for (i = 0; i < numCourseG; i++)
	{
		cout << "Course " << i + 1 << " is: " << courseListG[i].get_course_name() << endl;
		cout << "Course # is: " << courseListG[i].get_course_number() << endl << endl;
	}
}
