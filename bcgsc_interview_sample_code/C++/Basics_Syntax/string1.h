#ifndef _STRING_H_
#define _STRING_H_

class String
{

public:
  String();
  String(const String &s);
  String(const char *s);

  ~String();

  void  set_data(const String &s);
  void  set_data(const char *s);

  int   length() const { return dm_len; };

  const char*	get_data() const { return dm_data; };

private:
  void set_to_empty_string(); // utility method

  int		dm_len;
  char*	dm_data;
};

#endif
