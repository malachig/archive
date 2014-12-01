/***
*	This program is a test function for the course class.
*	Several course objects will be created as a linked list using the 
*	linked list class.  The usage of the linked list class is also illustrated.
*	The courses added to the list will be printed.  Then one will be deleted, 
*	and the resulting list printed again.  
*	Usage a static global linked list is also illustrated.
*	Note: Normally I would not include so much code in main that need not be there,
*	however this is only a preliminary test to prove that my Course class works.
***/

#include <iostream>
using namespace std;

#include "my_string.h"
#include "course.h"
#include "linked_list.h"

// Function Prototypes
void add_course_to_list(Course *c, LinkedList &list);
void print_list(const LinkedList &list);

// Global variable declaration.
static LinkedList gList;  // Creates an empty list.

// A function allowing any function to get a local reference to the master list.
LinkedList& get_master_list() {return gList;}

void main()
{
	char str[100];  // Temporary input storage.
	int i;			// Loop control variable.
	Course *c;		// Points to any course object.

	LinkedList &list = get_master_list();  // get local reference to master list.

	// Add 5 dynamic course objects as a simple test of the class.
	for (i=0; i<5; i++)
	{
		cout << "Enter a course name: ";
		cin.getline(str, 100);

		c = new Course(str);  // Create course using overloaded constructor.

		add_course_to_list(c, gList);  // passes Global list by name.
	}

	// Call function to get master list and pass on its return to print_list()
	print_list(get_master_list());

	/* Remove the 3rd entry from the list, note, object is NOT deleted.  It is 
	*  only removed from the list.  Note: the list keeps track of void* to data,
	*  therefore we must cast to the actual data type, pointers to Course in 
	*  this case ... */

	c = (Course*) list.remove(2);

	// Will not need it so now we can delete the course object.
	delete c;

	// Print the remaining courses in the list.
	print_list(list);

	// Memory Clean Up
	if (list.size() > 0)
	{
		cout << endl << list.size() << " Courses Still in Course List" << endl;
		cout << "Cleaning up Dynamic Memory" << endl;
		
		for (int i = 0; i < list.size(); i++)
		{
			c = (Course*)list.remove(i); // uses operator overloading for '[]'
			delete c;
		}
		cout << "Complete!" << endl;
	}
}


// Will be adding data to list, ie, will be Changing it, even if done via
// a method, accept a reference, NOT constant reference, so we can cahgne data.
void add_course_to_list(Course *c, LinkedList &list)
{
	if (c == 0)  // check for bad data
		return;

	list.add(c);  // keeps track of object.  Note that the linked list keeps track
				  // of it as a 'void' pointer, it will have to be cast to a 'Course'
				  // object when it is retrieved from the list at a later time.
}

// Accept a reference to a LinkedList to avoid making copies of the entire list
// don't plan to change it, so accept is as a const reference.
void print_list(const LinkedList &list)
{
	Course *c;  // Local pointer to course object

	cout << "\n\nPrint List\n\n";

	for (int i=0; i<list.size(); i++)
	{
		c = (Course*)list[i];  // uses operator overloading, the index '[]' operator.

		// Use overloaded '<<' operator to print Course class.  Note: 'c' is a course
		// pointer.  The operator has been defined to work on course objects.  Therefore
		// we must dereference the pointer, otherwise it would print the contents of the 
		// pointer, which is an address.

		cout << i << ". - " << (*c);

		cout << endl;
	}
}

