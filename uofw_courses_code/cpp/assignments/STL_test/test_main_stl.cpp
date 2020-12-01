#include <iostream>
#include <list>
using namespace std;

#include "my_string.h"
#include "person.h"

// Function Prototypes
void add_person_to_list( Person *p, list<Person*> &pList );
// note, can't make list const !!
void print_list( list<Person*> &pList );
list<Person*>& master_person_list();
list<Person*>::iterator find_person(list<Person*> &pList, Person *pPtr);

static list<Person*>  pList_g; // creates the default empty list, ie, size is zero


int main()
{
  char str1[100];
  char str2[100];
  string name;   // for input
  string addr;
  int     i;
  list<Person*>::iterator test;  // simple interator variable for a person.
  list<Person*>::iterator personFound;
  Person  *p;    // points to any Person object
  Person *p2;
  Person *p3;
  int size;

  list<Person*> &pList = master_person_list();

  // add 3 *dynamic* objects
  for ( i=0; i<3; i++ ) 
  {
    cout << "Enter a name: ";
    cin.getline(str1,100);
	name = str1;

	cout << "Enter an Address: ";
    cin.getline(str2,100);
	addr = str2;
  
    p = new Person(str1, str2);  // create dynamic object
  
    add_person_to_list( p, pList );
  }

  p2 = new Person("dummy1", "fuckville");
  add_person_to_list(p2, pList);
  p3 = new Person("dummy2", "shithole");
  add_person_to_list(p3, pList);

  print_list(pList_g);

  // Determining size of the list
  size = pList.size();
  cout << "\nThe size of the list is: " << size << endl;


  // Try removing an object from the middle of the list
	personFound = find_person(pList, p2);
	test = pList_g.erase(personFound);

	cout << "\nAfter removing an element from inside the list:" << endl;
	print_list(pList);


  // Memory Clean-up Attempt
  list<Person*>::iterator itr = pList.begin();
  
  while (itr != pList.end())
  {
	p = *itr;
	delete p;
	itr++;
  }
    
  // clear the entire list, note, it does NOT delete the data, so we now have a memory
  // leak because we lost our only connection to the dynamic objects
  pList.clear();

  print_list( pList ); 

  return 0;
}


// will be adding data to list, ie, will be CHANGING it, even if done via
// a method, accept a reference, NOT constant reference, so we can change data
void add_person_to_list(Person *pPtr, list<Person*> &pList )
{
  if ( pPtr == 0 )  // check for bad data
    return;

  pList.push_back( pPtr ); // keeps track of Person* data
}



// accept a refernce to LinkedList to avoid making copies of the entire list
// don't plan to change it, so SHOULD accept it as a const reference
// HOWEVER the iterator is non-const pecasue you can do '++' to it
// sooo, list *has* to be non-const... C++ at work for ya
void print_list( list<Person*> &pList )
{
  const Person* pPtr;  // local pointer to String object
  list<Person*> &local_list = master_person_list();

  cout << "\n\n Printing List\n\n";

  if (pList.empty())
  {
	  cout << "\nThe list is empty" << endl;
	  return;
  }

  // use iterator to 'loop' over list
  list<Person*>::iterator itr = pList.begin();

  while ( itr != pList.end() ) 
  {
    // list keeps track of String* data, no casting needed
    // uses operator overlaoding, '*' actually returns an element in list
    pPtr = *itr; 

    cout << *itr << endl;
	pPtr->print();

    itr++;  // advance to next element in list
  }

}

list<Person*>& master_person_list()
{
	return pList_g;
}

list<Person*>::iterator find_person(list<Person*> &pList, Person *pPtr)
{
  Person  *p;
  list<Person*>::iterator itr1 = pList.begin();
  
  while (itr1 != pList.end())
  {
	  p = *itr1;
	  if (p == pPtr)
		return itr1;

	  itr1++;
  }
  return 0;
}