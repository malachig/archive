#include <iostream>
using namespace std;

// default values specified in the prototype
void foo (short p1, short p2=0, short p3=0);

// default values *not* specified in function code
void foo (short a, short b, short c)
{
	cout << "foo(): " << a << '\t' << b << '\t' << c << endl;
}

int main()
{
	foo(1);		// output -> foo(): 1 0 0
	foo(1,2);	// output -> foo(): 1 2 0
	foo(1,2,3);	// output -> foo(): 1 2 3

	return 0;
}
