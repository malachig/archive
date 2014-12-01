/***
*	The method code (definition) for the Student Class.
*	Located in a seperate .cpp file for clarity and modularity. 
*	The class decleration (its name, methods and data fields) are found
*	in the student.h header file and included when needed.
***/

// So compiler knows class declaration (methods and data members to expect)
#include "student.h"	

// Note: Method names must be scoped to the class they belong.

/***  STUDENT CLASS METHOD CODE ***/

// Default Constructor.
Student::Student()
{
	dm_studentNumber = 0;
}


// Copy Constructor
Student::Student(const Student &s)
{
	dm_studentNumber = s.dm_studentNumber;
}


// Overloaded Constructor.
// This one allows user to initialize Student attributes when the student
// object is created.
Student::Student(unsigned long number)
{
	dm_studentNumber = number;
}


Student::~Student()
{
	// Does nothing for now, but defined just for kicks.
}


// Overloading the assignment '=' operator.
const Student& Student::operator=(const Student &s)
{
	if (this == &s)		// If LHS is the same as RHS
		return *this;	
	
	dm_studentNumber = s.dm_studentNumber;
	dm_courseList = s.dm_courseList;

	return *this;
}


// Method Code for Accessors

unsigned long 
Student::get_student_number() const
{
	return dm_studentNumber;
}


void 
Student::print() const
{
	Course *c;				// Local pointer to course object.

	// Print out the person's name and address.
	Person::print();

	// Print out the student information.
	cout << "STUDENT INFO" << endl;
	// Now actually display the info to the screen.
	cout << "The student # is: " << dm_studentNumber << endl;

	if (dm_courseList.size() < 1)
	{
		cout << "\nStudent has no courses\n" << endl;
		return;
	}

	// Display list of courses for the student.
	cout << "\nThe student is enrolled in the following courses:" << endl;
	
	for (int i = 0; i < dm_courseList.size(); i++)
	{
		c = (Course*)dm_courseList[i];  // Get a course from the list.

		cout << i+1 << ". - " << (*c);
		cout << endl;
	}
	cout << endl;
}


void 
Student::print_courses () const
{
	Course *c;  // Local pointer to a course object.

	cout << "That student is enrolled in the following courses:" << endl;
	
	for (int i = 0; i < dm_courseList.size(); i++)
	{
		c = (Course*)dm_courseList[i];  // Get a course from the list.

		cout << i+1 << ". - " << (*c);
		cout << endl;
	}
}


void 
Student::add_course(Course *c)
{
	// Check for bad data.
	if (c == 0)  
		return;

	// Check to see if that course is already in the student's list.
	if (dm_courseList.size() > 0)  // no point if the list is empty.
	{
		// Use Object Address Comparison to see if the selected course is 
		// already in the student's list.
		for (int i = 0; i < dm_courseList.size(); i++)
		{
			if (dm_courseList[i] == c)
			{
				cerr << "\n*** Student already enrolled in that course ***" << endl;
				return;
			}
		}
	}

	dm_courseList.add(c);  // Add the course to the list.
	cout << "Course Added" << endl;
}


void 
Student::remove_course(Course *c)
{
	Course *removeCourse;
	
	if (c == 0)  // Check for bad data.
		return;

	// Check if the course to be removed is in the list.
	if (dm_courseList.size() > 0)  // Make sure the list is not empty.
	{
		for (int i = 0; i < dm_courseList.size(); i++)
		{
			// Use Object Address Comparison to remove the target course.
			if (dm_courseList[i] == c)
			{
				// Remove course from list and Update the linked list 
				removeCourse = (Course*)dm_courseList.remove(i);

				// Do not delete the memory associated with this pointer now because
				// the student does not OWN the course objects.  It has a list of 
				// nodes which point to courses in the master course list.  Therefore
				// numerous other students could be using those same course objects.
				
				// Exit function once the target course has been found and removed.
				return;
			}
		}
		// If the course was not found display an error message
		cerr << "\n*** Student not enrolled in that course ***" << endl;
	}
}

/***  GRAD_STUDENT CLASS METHOD CODE ***/

// Default Constructor.
GradStudent::GradStudent()
{
	dm_thesisTitle = "Thesis Title Unknown";
	dm_advisorName = "Advisor Name Unknown";
}


// Copy Constructor
GradStudent::GradStudent(const GradStudent &gs)
:dm_thesisTitle(gs.dm_thesisTitle), dm_advisorName(gs.dm_advisorName)
{
	
}

// Overloaded Constructor.
// This one allows user to initialize Grad Student attributes when the 
// object is created.
GradStudent::GradStudent(const String &thesis, const String &advisor)
: dm_thesisTitle(thesis), dm_advisorName(advisor)
{

}


GradStudent::~GradStudent()
{
	// Does nothing for now, but defined just for kicks.
}


// Overloading the assignment '=' operator.
const GradStudent& GradStudent::operator=(const GradStudent &gs)
{
	if (this == &gs)		// If LHS is the same as RHS
		return *this;	
	
	dm_thesisTitle = gs.dm_thesisTitle;
	dm_advisorName = gs.dm_advisorName;
	
	return *this;
}


// Method Code for Mutators
void 
GradStudent::set_thesis_title(const String &thesis)
{
	dm_thesisTitle =thesis;
}


void 
GradStudent::set_advisor_name(const String &advisor)
{
	dm_advisorName = advisor;
}


void 
GradStudent::print() const
{
	// Print out student and person info
	Student::print();

	// Now actually display the Grad info to the screen.
	cout << "GRAD INFO" << endl;
	cout << "The Grad Student's Thesis Title is:\n" << dm_thesisTitle << endl;
	cout << "The Name of their Advisor is: " << dm_advisorName << endl;
	cout << endl;
}


/***  UNDER_GRAD_STUDENT CLASS METHOD CODE ***/

// Default Constructor.
UnderGradStudent::UnderGradStudent()
{
	dm_projectTitle = "Project Title Unknown";
}


// Copy Constructor
UnderGradStudent::UnderGradStudent(const UnderGradStudent &ugs)
:dm_projectTitle(ugs.dm_projectTitle)
{
	
}

// Overloaded Constructor.
// This one allows user to initialize Under Grad Student attributes when the 
// object is created.
UnderGradStudent::UnderGradStudent(const String &project)
: dm_projectTitle(project)
{

}


UnderGradStudent::~UnderGradStudent()
{
	// Does nothing for now, but defined just for kicks.
}


// Overloading the assignment '=' operator.
const UnderGradStudent& UnderGradStudent::operator=(const UnderGradStudent &ugs)
{
	if (this == &ugs)		// If LHS is the same as RHS
		return *this;	
	
	dm_projectTitle = ugs.dm_projectTitle;
		
	return *this;
}


// Method Code for Mutators
void 
UnderGradStudent::set_project_title(const String &project)
{
	dm_projectTitle = project;
}


void 
UnderGradStudent::print() const
{
	// Print out student and person info
	Student::print();

	// Now actually display the UnderGrad info to the screen.
	cout << "UNDERGRAD INFO" << endl;
	cout << "\nThe UnderGrad Student's Project Title is:\n" << dm_projectTitle << endl;
	cout << endl;
}
