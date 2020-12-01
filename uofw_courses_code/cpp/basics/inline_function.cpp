#include <iostream>
using namespace std;

inline int maxValue (int a, int b);

int main()
{
	int a, b;

	cin >> a >> b;
	cout << maxValue(a,b) << endl;

	return 0;
}

int maxValue (int a, int b)
{
	if (a > b)
		return a;
	return b;
}
