#include <iostream>
#include <cstring>
#include <cstdlib>
using namespace std;

#define SIZE 25

class strtype
{
private:
	char *p;
	int len;

public:
	strtype();	// constructor
	~strtype();	// destructor
	void set(char *ptr);
	void show();
};

// Initialize a string object.
strtype::strtype()
{
	p = (char *) malloc (SIZE);

	if(!p)
	{
		cout << "Allocation error\n";
		exit(1);
	}

	*p = '\0';
	len = 0;
}

// Freee memory when destroying string object.
strtype::~strtype()
{
	cout << "Freeing p\n";
	free(p);
}

void strtype::set(char *ptr)
{
	if(strlen(p) >= SIZE)
	{
		cout << "String too big\n";
		return;
	}
	strcpy(p, ptr);
	len = strlen(p);
}

void strtype::show()
{
	cout << p << " - length: " << len;
	cout << "\n";
}

int main()
{
	strtype s1, s2;

	s1.set("this is a test.");
	s2.set("I like C++.");

	s1.show();
	s2.show();

	return 0;
}


		

