// Header
// Title: Comp finance - Mini task 2
// Student ID: 10134621
// Date Created: 03/03/21
// Last Edited: 03/03/21


#define _USE_MATH_DEFINES_


// Includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <fstream>
#include <vector>
#include <constants.h>  // header file for constants


// Decalre Functions

// calculate cummulative normal distribution
double norm_cum(const double& x);

// calculate f
double f(const double& r, const double& t, const double& T);

// calculate m
double m(const double& r, const double& t, const double& T);

// calculate q
double q(const double& t, const double& T);

// calculate v^2
double v2(const double& t, const double& T);

// calculate P
double P(const double& r, const double& t, const double& T);

// calculate n
double n(const double& r, const double& t, const double& T);

// calculate k
double k2(const double& t, const double& T);

// calculate V for put
double V_put(const double& r, const double& t, const double& T, const double& h);



// Begin main program
int main() 
{
	// define variables
	const double t{ 0 };
	const double T{ 3 };
	double b = 0.2;  // lower r limit
	double a = 0;  // upper r limit
	double number_calc = 100;  // number of calculations

	// open a file stream for writing
	std::ofstream output;

	// open the csv file
	output.open("data.csv");

	// if the file is open
	if (output.is_open()) {
	
		// for loop over r values
		for (double r{ 0 }; r <= b + 0.002; r += (b - a) / number_calc) {

			// calculate h
			double h_val = (constants::X_r - f(r, t, T)) / pow(v2(t, T), 0.5);

			// calculate P
			double P_val = P(r, t, T);

			// calculate V(r, t=0, T) for a put option
			double V_val = V_put(r, t, T, h_val);

			// write data to file
			output << r << "," << P_val << "," << V_val << std::endl;

			// print data to screen
			std::cout << "r = " << r << ", P = " << P_val << ", V_val = " << V_val << std::endl;

		}

		// close the file
		std::cout << "File write successful" << std::endl;
		output.close();

	}
	// if file could not be opened
	else {
		std::cout << "Error: could not open file" << std::endl;
		return 1;
	}

	// Test r_0 values
	std::cout << m(constants::r_0, t, T) << std::endl;
	std::cout << n(constants::r_0, t, T) << std::endl;
	std::cout << k2(t, T) << std::endl;
	std::cout << q(t, T) << std::endl;
	std::cout << f(constants::r_0, t, T) << std::endl;
	std::cout << v2(t, T) << std::endl;
	std::cout << (constants::X_r - f(constants::r_0, t, T)) / pow(v2(t, T), 0.5) << std::endl;
	std::cout << norm_cum((constants::X_r - f(constants::r_0, t, T)) / pow(v2(t, T), 0.5)) << std::endl;
	std::cout << P(constants::r_0, t, T) << std::endl;
	std::cout << V_put(constants::r_0, t, T, (constants::X_r - f(constants::r_0, t, T)) / pow(v2(t, T), 0.5)) << std::endl;


	return 0;
}  // End main program


// Define functions

// calculate V for put
double V_put(const double& r, const double& t, const double& T, const double& h) 
{
	return P(r, t, T) * norm_cum(h);
}

// calculate cummulative normal distribution
double norm_cum(const double& x) 
{
	return 0.5 * erfc(-x / pow(2, 0.5));
}

// calculate f
double f(const double& r, const double& t, const double& T) 
{
	return m(r, t, T) - 0.5 * q(t, T);
}

// calculate m
double m(const double& r, const double& t, const double& T) 
{
	return exp(-constants::kappa*(T-t))*r+(1-exp(-constants::kappa*(T-t)))*constants::theta;
}

// calculate q
double q(const double& t, const double& T) 
{
	return (pow(constants::sigma, 2) / (3 * pow(constants::kappa, 2)))* pow(1 - exp(-constants::kappa * (T - t)), 5);
}

// calculate v^2
double v2(const double& t, const double& T) 
{
	return (pow(constants::sigma, 2) / constants::kappa) * (1 - exp(-constants::kappa * (T - t)));
}

// calculate P
double P(const double& r, const double& t, const double& T) 
{
	return exp((2. / 3.) * k2(t, T) - (1. / 4.) * n(r, t, T));
}

// calculate n
double n(const double& r, const double& t, const double& T) 
{
	return r * (T - t) - ((constants::theta - r) / (2 * constants::kappa)) * (1 - exp(-4 * constants::kappa * (T - t)));
}

// calculate k
double k2(const double& t, const double& T) 
{
	return ((pow(constants::sigma, 2)) / (2 * pow(constants::kappa, 3))) * (5 * exp(-constants::kappa * (T - t)) - 3 * exp(-2 * constants::kappa * (T - t))
		+ 3 * constants::kappa * (T - t) - 2);
}