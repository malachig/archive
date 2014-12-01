#include <iostream>
using namespace std;

int main()
{
	char *name = "Frodo";

	cout << "Char: " << name[0] << endl;
	cout << "Short: " << (short)(name[0]) << endl;
	cout << "String: " << name << endl;
	cout << "Address: " << (unsigned long) name << endl;

	return 0;
}
