#include <iostream>
using namespace std;

#include "my_string.h"

#define errMsg "<Error> String::"

ostream& operator<<(ostream &os, const String &s)
{
  os << s.dm_data;
  return os;
}



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
  If pointer is null, creates the empty string.
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
	Assignment operator.
*****/
const String& String::operator=(const String &s)
{
  if ( this == &s )	// same object, do not copy
    return *this;

  set_data(s.dm_data);

  return *this;
}

/****
  Assignment operator. Does not create temporary String if argument is char array.
*****/
const String& String::operator=(const char *s)
{
  set_data(s);

  return *this;
}


/****
  set_data() - makes deep copy of given character array.
  Leaves string as is if pointer is null.
*****/
void String::set_data(const String &s)
{
  set_data( s.dm_data );
}

/****
  set_data() - makes deep copy of given character array.
  Leaves string as is if pointer is null.
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
char& String::operator[]( int idx )
{
  if ( idx<0 || idx>=dm_len )
  {
    static char tmp;	// static so it doesn't go away...
    cerr << errMsg << "operator[]" << endl;
    cerr << " out of bounds:"<<idx<<endl;
    return tmp; 
  }

  return dm_data[idx];
}


/****

*****/
char String::operator[]( int idx ) const
{
  if ( idx<0 || idx>=dm_len )
  {
    cerr << errMsg << "operator[]" << endl;
    cerr << " out of bounds:"<<idx<<endl;
    return '\0';
  }

  return dm_data[idx];
}


/****

*****/
void String::set_to_empty_string()
{
  dm_len = 0;
  dm_data = new char[1];
  dm_data[0] = '\0';
}



