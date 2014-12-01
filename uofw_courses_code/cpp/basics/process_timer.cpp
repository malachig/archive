#include <iostream>
#include <ctime>
using namespace std;

class timer 
{
private:
	clock_t start;

public:
	timer();	// constructor
	~timer();	// destructor
};

timer::timer()
{
	start = clock();
}

timer::~timer()
{
	clock_t end;

	end = clock();

	cout << "Elapsed time: " << (end - start) / CLOCKS_PER_SEC << endl;
}


int main()
{
	timer obj;
	char c;

	// delay ...

	cout << "Press a key followed by ENTER: ";
	cin >> c;

	return 0;
}
