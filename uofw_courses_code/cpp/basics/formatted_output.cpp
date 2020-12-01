#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

int main()
{
	float amount, principal = 1000.0, rate = 0.05;

	cout << "Year" << setw(21) << "Amount on deposit" << endl;

	for (int yr = 1; yr <= 10; yr++)
	{
		// a = p(1+r)n
		amount = principal * pow(1+rate, yr);

		cout << setw(4) << yr << setiosflags(ios::fixed|ios::showpoint)
			 << setw(21) << setprecision(2) << amount << endl;
	}

	return 0;
}
