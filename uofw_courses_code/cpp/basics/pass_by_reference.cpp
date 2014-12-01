#include <iostream>
using namespace std;

//function prototype using inline method to place the code right away
void double_value (int &val)
{
	val *= 2;
}

int main()
{
	int num;

	double_value(num);

	return 0;
}