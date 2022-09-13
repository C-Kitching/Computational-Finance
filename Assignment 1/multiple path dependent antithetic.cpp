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
	double K{ 35 };  // points in the sample path

	// containers
	std::vector<double> N_store;
	std::vector<std::vector<double>> master_time_store;

	// loop over number of calculations
	for (int i{ 0 }; i < 100; i++) {

		std::vector<double> time_store;
		
		for (int N{ 100000 }; N <= 1500000; N += 100000) {

			if (i == 0) N_store.push_back(N);

			auto start = std::chrono::steady_clock::now();  // get start time
			double value = value_Asian_call(initial_share_price, interest_rate, dividend_rate, volatility, expiration, N, K);
			auto finish = std::chrono::steady_clock::now();  // get finish time
			auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>> (finish - start);  // convert into seconds

			time_store.push_back(elapsed.count());
		}

		master_time_store.push_back(time_store);
	}

	// calculate average
	std::vector<double> average;
	for (int i{ 0 }; i < N_store.size(); i++) {

		double sum{ 0 };
	
		for (int j{ 0 }; j < master_time_store.size(); j++) {
			sum += master_time_store[j][i];
		}

		average.push_back(sum / master_time_store.size());	
	}

	// open a file stream for writing
	std::ofstream output;

	// open the csv file
	output.open("path dep antithetic.csv");

	// if the file is open
	if (output.is_open()) {

		// loop over data containers
		for (int i{ 0 }; i < N_store.size(); i++) {

			// write data to file
			output << N_store[i] << "," << average[i] << std::endl;
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

		// time step
		double dt{ expiration / K };

		// containers for stock path
		std::vector<double> stock_path1;
		std::vector<double> stock_path2;

		// record initial values
		stock_path1.push_back(initial_share_price);
		stock_path2.push_back(initial_share_price);

		// generate stock path
		for (int i{ 1 }; i <= K; i++) {

			// generate random number
			double phi = ND(rnd);

			// gemerate stock path
			stock_path1.push_back(stock_path1[i - 1] * exp((interest_rate - dividend_rate -
				0.5 * pow(volatility, 2)) * dt + volatility * phi * pow(dt, 0.5)));
			stock_path2.push_back(stock_path2[i - 1] * exp((interest_rate - dividend_rate -
				0.5 * pow(volatility, 2)) * dt - volatility * phi * pow(dt, 0.5)));
		}

		// calculate A
		double A1, A2;
		double A1_sum{ 0 }, A2_sum{ 0 };
		for (int i{ 1 }; i <= K; i++) {
			A1_sum += stock_path1[i];
			A2_sum += stock_path2[i];
		}
		A1 = A1_sum / K;
		A2 = A2_sum / K;

		// add in the payoff
		sum += std::max(stock_path1.back() - A1, 0.);
		sum += std::max(stock_path2.back() - A2, 0.);
	}

	// average over all paths
	return exp(-interest_rate * expiration) * sum / (2.*N);
}