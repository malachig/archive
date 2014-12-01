#include <iostream>
using namespace std;

#include "string1.h" // Have changed name because my Visual C++ Studio
					 // already has a string.h and it was being included
					 // by default.


#define errMsg "<Error> String::" 


/****
  Default constructor. Creates an empty string.
*****/
String::String()
{
  set_to_empty_string();
}


/****
  Copy constructor. Creates a deep copy of the given string.
*****/
String::String(const String &s)
{
  if ( s.dm_data != NULL ) {
    dm_len = s.dm_len;
    dm_data = new char[dm_len+1];
    strcpy(dm_data, s.dm_data);
  }
  else {
    set_to_empty_string();
  }
}



/****
  Constructor. Makes a deep copy of given character array.
  If pointer argument is null, creates the empty string.
*****/
String::String(const char* s)
{
  if ( s == NULL ) 
  {
    cerr << errMsg << "String(const char*)" << endl;
    cerr << " NULL Pointer."<<endl;

    set_to_empty_string();
  }
  else
  {
    dm_len = strlen(s);
    dm_data = new char[dm_len+1];
    strcpy(dm_data, s);
  }
}


/****
  Destructor.
*****/
String::~String()
{
  delete [] dm_data;
}


/****
  set_data() - makes deep copy of given character array.
  Leaves string as is if pointer is null.
*****/
void String::set_data(const String &s)
{
  set_data(s.dm_data);
}


/****
  set_data() - makes deep copy of given character array.
  Leaves string data member unchanged if pointer argument is null.
*****/
void String::set_data(const char *s)
{
  if ( s == NULL ) {
    cerr << errMsg << "set_data(const char *s)" << endl;
    cerr << " NULL Pointer."<<endl;
  }
  else
  {
    dm_len = strlen(s);
    delete [] dm_data;
    dm_data = new char[dm_len+1];
    strcpy(dm_data, s);
  }
}



/****
  
*****/
void String::set_to_empty_string()
{
  dm_len = 0;
  dm_data = new char[1];
  dm_data[0] = '\0';
}



