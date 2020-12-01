#include <iostream>
using namespace std;

/***
* When instructions are identical and functions only differ in data type,
* you can use template functions instead of function overloading.
***/

template <class T>
void swap( T &a, T &b)
{
	T tmp;

	tmp = a;
	a = b;
	b = tmp;
}

template <class T>
void bubble_sort(T list[], int size)
{
	int i;
	int flg = TRUE;

	while (flg = TRUE)
	{
		flg = FALSE;

		for (i = 0; i < size-1; i++)
		{
			if(list[i+1] < list[i])
			{
				swap(list[i], list[i+1]);
				flg = TRUE;
			}
		}
	}
}

/***
* When bubble_sort is called like this:
*	int arr[100];
*	bubble_sort(arr, 100);
* The compiler will generate the code for bubble_sort replacing every 
* occurance of T with int.  Same thing will happen with swap.
***/