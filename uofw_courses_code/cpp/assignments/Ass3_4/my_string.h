#ifndef __MY_STRING_HPP
#define __MY_STRING_HPP

#include <iostream>
using namespace std;

class String; // forward class decleration

// M$ C++ wants this function prtotype declared *outside* of friend class as well but
// an argument is a reference to a String, thus need the previous forward class decleration
ostream& operator<<(ostream &os, const String &s);


class String
{
  friend ostream& operator<<(ostream &os, const String &s);

public:
  String();
  String(const String &s);
  String(const char *s);
  ~String();

  const String& operator=(const String &s);
  const String& operator=(const char *s);

  char& operator[]( int idx );
  char  operator[]( int idx ) const;

  int length() const { return dm_len; };

  const char*	data() const { return dm_data; };

  void  set_data(const char *s);
  void  set_data(const String &s);

private:
  void set_to_empty_string();

  int		dm_len;
  char*	dm_data;
};


#endif
