#include <iostream>
using namespace std;

class myclass {
	int a;
public:
	myclass();	// constructor
	~myclass(); // destructor
	void show();
};

myclass::myclass()
{
	cout << "In constructor\n";
	a = 10;
}

myclass::~myclass()
{
	cout << "Destructing...\n";
}

void myclass::show()
{
	cout << a << "\n";
}



int main()
{
	myclass obj;

	obj.show();

	return 0;
}
