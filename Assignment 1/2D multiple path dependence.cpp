// HEADER
// Student ID: 10134521
// Title: Assignment 1
// Date Created: 18/03/21
// Last Edited:


// math constants
#define _USE_MATH_DEFINES_


// Includes
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <vector>
#include <chrono>


// Function declerations

// value Asian call
double value_Asian_call(const double& initial_share_price, const double& interest_rate, const double& dividend_rate, const double& volatility,
	const double& expiration, const int& N, const double& K);

// Begin main program
int main()
{
	// define parameters
	double expiration{ 1.25 };
	double volatility{ 0.37 };
	double interest_rate{ 0.03 };
	double dividend_rate{ 0.04 };
	double initial_share_price{ 900 };

	// containers to store data
	std::vector<int> parameter_store;
	std::vector<double> value_store;
	std::vector<double> time_store;

	int N{ 200000 };

	// loop over different K
	for (int K{ 1 }; K <= 50; K += 1) {

		// store N and K
		parameter_store.push_back(K);

		auto start = std::chrono::steady_clock::now();  // get start time

		// calculate and store option value 
		value_store.push_back(value_Asian_call(initial_share_price, interest_rate, dividend_rate, volatility, expiration, N, K));

		auto finish = std::chrono::steady_clock::now();  // get finish time
		auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>> (finish - start);  // convert into seconds
		time_store.push_back(elapsed.count());  // store the time
	}

	// open a file stream for writing
	std::ofstream output;

	// open the csv file
	output.open("path dependence.csv");

	// if the file is open
	if (output.is_open()) {

		// loop over data containers
		for (int i{ 0 }; i < parameter_store.size(); i++) {

			// write data to file
			output << parameter_store[i] << "," << value_store[i] << "," << time_store[i] << std::endl;
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


	return 0;
}  // End main progrma


// Function definitions

// value Asian call
double value_Asian_call(const double& initial_share_price, const double& interest_rate, const double& dividend_rate, const double& volatility,
	const double& expiration, const int& N, const double& K)
{
	// declare random number generator
	static std::mt19937 rnd;

	// declare the normal distrubtion
	std::normal_distribution<double> ND(0., 1.);

	// initalise sum to zero
	double sum{ 0 };

	// loop over all MC paths
	for (int i{ 0 }; i < N; i++) {

		// create a sample path
		double dt{ expiration / K };
		std::vector<double> stock_path;
		stock_path.push_back(initial_share_price);
		for (int i{ 1 }; i <= K; i++) {

			// generate random number
			double phi = ND(rnd);

			// gemerate stock path
			stock_path.push_back(stock_path[i - 1] * exp((interest_rate - dividend_rate -
				0.5 * pow(volatility, 2)) * dt + volatility * phi * pow(dt, 0.5)));
		}

		// calculate A
		double A;
		double A_sum{ 0 };
		for (int i{ 1 }; i <= K; i++) A_sum += stock_path[i];
		A = A_sum / K;

		// add in the payoff
		sum += std::max(stock_path.back() - A, 0.);
	}

	// average over all paths
	return exp(-interest_rate * expiration) * sum / N;
}