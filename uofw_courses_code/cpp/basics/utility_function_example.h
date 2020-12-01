#ifndef _UTILITY_FUNCTION_EXAMPLE.H   // to avoid multiple and recursive inclusions
#define _UTILITY_FUNCTION_EXAMPLE.H

class SalesPerson
{
public:
	SalesPerson();

	void set_monthly_sales(int month, float amount);

	void print_annual_sales() const;

private:
	float total_annual_sales() const;	// utility function

	// ** data members ** //
	float dm_sales[12];		// total sales for each month
};

// ...

void SalesPerson::print_annual_sales() const
{
	// ** may loop over months & print monthly sales ** //

	cout << setprecision(2)
		<< setiosflags(ios::fixed||ios::showpoint)
		<< endl << "The total annual sales are: $"
		<< total_annual_sales() << endl;  // calls private utlility function total_annual_sales().
}
