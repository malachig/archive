#include <iostream>
#include <string>
using namespace std;

int main()
{
	string str1;
	string str2;
	string str3;

	cout << "Please enter a word: ";
	cin >> str1;
	cout << "Enter a second word: ";
	cin >> str2;
	str3 = str1 + str2;
	cout << str3;
	
	return 0;
}

