#include <iostream>
using namespace std;

class myclass
{
private:
	int a;

public:
	myclass(int x);	// constructor
	void show();
};

myclass::myclass(int x)
{
	cout << "In constructor\n";
	a = x;
}

void myclass::show()
{
	cout << a << "\n";
}

int main()
{
	myclass obj(4);

	obj.show();

	return 0;
}
