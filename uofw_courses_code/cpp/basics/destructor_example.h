#ifndef _DESTRUCTOR_EXAMPLE_H_
#define _DESTRUCTOR_EXAMPLE_H_

class IntVec
{
public:
	IntVec();
	IntVec(const IntVec &vec);
	IntVec(int lengt, int val);
	~IntVec;  // destructor.

	void resize(int newLength);
	int length() const {return dm_length};
	void set_element(int idx, int val);
	int element(int idx) const;

private:
	int* dm_vec;	// variable length vector.
	int dm_length;
};

// *** Implementation *** //

IntVec::IntVec()
{
	dm_length = 0;
	dm_vec = null;
}

// copy constructor
IntVec::IntVec(const IntVec &vec)
{
	dm_vec = null;
	dm_length = vec.dm_length;

	if (dm_length > 0)
	{
		dm_vec = new int[dm_length];
		for (int i = 0; i < dm_length; i++)
			dm_vec[i] = vec.dm_vec[i];
	}
}

// Overloaded constructor
IntVec::IntVec(int length, int val)
{
	if ((dm_vec = new int[length]) != null)
		dm_length = length;
	else
		dm_length = 0;

	for (int i = 0; i < dm_length; i++)
		dm_vrc[i] = val;
}

IntVec::~IntVec()
{
	delete [] dm_vec;  // if null, it's safe.
					   // members dm_length & dm_vec not destroyed yet
}

void IntVec::resize(int newLength)
{
	if (newLength <= 0)
	{
		delete [] dm_vec;
		dm_length = 0;
		return;
	}

	int *newVec = new int[newLength];
	int ln = dm_length;

	if (newLength < ln)
		ln = newLength;

	for (int i=0; i<ln; i++)
		newVec[i] = dm_vec[i];

	delete [] dm_vec;
	dm_vec = newVec;
}

void IntVec::set_element(int idx, int val)
{
	if (idx >= 0 && idx < dm_length)
		dm_vec[idx] = val;
}

int IntVec::element(int idx) const
{
	int val = 0;  // what should be returned if out of bounds?
				  // if anything? exception handling (later)

	if (idx >= 0 && idx < dm_length)
		val = dm_vec[idx];

	return val;
}


























#endif