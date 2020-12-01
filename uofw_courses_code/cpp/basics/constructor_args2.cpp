#include <iostream>
using namespace std;

class myclass
{
private:
	int a, b;

public:
	myclass(int x, int y);	// constructor
	void show();
};

myclass::myclass(int x, int y)
{
	cout << "In constructor\n";
	a = x;
	b = y;
}

void myclass::show()
{
	cout << a << ' ' << b << "\n";
}

int main()
{
	myclass obj(4, 7);

	obj.show();

	return 0;
}

