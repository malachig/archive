#include <iostream>
using namespace std;

#include "linked_list.h"
#include "my_string.h"

void add_string_to_list( String *s, LinkedList &list );

void print_list( const LinkedList &list );


static LinkedList  gList; // creates the default empty list, ie, size is zero


LinkedList& get_master_list() { return gList; }


void main()
{
  char    str[100];   // for input
  int     i;
  String  *s;    // points to any String object

  LinkedList &list = get_master_list(); // get local reference to master list

  // add 5 *dynamic* objects
  for ( i=0; i<5; i++ ) 
  {
    cout << "Enter a string: ";
    cin.getline(str,100);
  
    s = new String( str );  // create dynamic object
  
    add_string_to_list( s, gList );	// passes global list by name
  }

  // calls function to get master list, passes on what it returns
  print_list( get_master_list() ); 	

  // remove the 3rd entry from list, note, object NOT deleted, just removed from the 
  // list, the list keeps track of void* to data, must cast to actual data type,
  // pointers to String in this case... 
  s = (String*) list.remove(2); 

  // will not need it, so now safe to delete the String object
  delete s;

  print_list( list );  // uses local reference to master list

  // clear the entire list, note, it does NOT delete the data, so we now have a memory
  // leak because we lost our only connection to the dynamic objects
  list.clear_list();

  print_list( list ); 
}

// will be adding data to list, ie, will be CHANGING it, even if done via
// a method, accept a reference, NOT constant reference, so we can change data
void add_string_to_list( String *sPtr, LinkedList &list )
{
  if ( sPtr == 0 )  // check for bad data
    return;

  list.add( sPtr ); // keeps track of object, note, linked list keeps track of it
                    // as 'void' pointer, will have to be cast to a 'String' object
                    // when it is retrieved from the list at a later time
}



// accept a refernce to LinkedList to avoid making copies of the entire list
// don't plan to change it, so accept it as a const reference
void print_list( const LinkedList &list )
{
  String *sPtr;  // local pointer to String object

  cout << "\n\n Printing List\n\n";

  for ( int i=0; i<list.size(); i++ ) 
  {
    // list keeps track of void* to data, must cast to actual data type,
    // pointers to String in this case... very dangerous way of doing things
    // but until learn about templates, the most convenient way to go.

    sPtr = (String*)list[i]; // uses operator overlaoding, the index '[]' operator 

    // use overloaded '<<' opertor to print String class, note, sPtr is a String pointer
    // the operator has been defined to work on String *objects*, so must dereference the
    // pointer, otherwise it would print the contents of the pointer, which is an address
    cout << i << ". - " << (*sPtr); 
    
    cout << endl;
  }

}