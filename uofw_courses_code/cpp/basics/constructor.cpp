#include <iostream>
using namespace std;


class myclass {
private:
	int i, j;

public:
	myclass(int a, int b);
	void show();
};


myclass::myclass(int a, int b)
{
	i = a;
	j = b;
}


void myclass::show()
{
	cout << i << ' ' << j << "\n";
}


int main()
{
	int x, y;

	cout << "Enter two integers: ";
	cin >> x >> y;

	// Use variables to construct object obj.
	myclass obj(x, y);

	obj.show();

	return 0;
}

