#include <iostream>
using namespace std;

// use const variables over #define since type checking 
// can be done.
// const variables can be initialized when defined.

const short kMaxNameLength = 40;

int main()
{
	char name[kMaxNameLength];
	short aShort;
	long aLong;
	float aFloat;

	cout << "Enter your name: ";

	cin >> name;
		// scans stream until after user presses <return>
		// until a white space is encountered
		// the rest is left on the stream
	cout << "Enter a short, long, and a float: ";
	cin >> aShort >> aLong >> aFloat;

	return 0;
}
