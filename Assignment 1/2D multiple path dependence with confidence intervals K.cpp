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
	std::vector<double> upper_confidence;
	std::vector<double> lower_confidence;

	// number of calculations
	int M{ 100 };

	// number of monte carlo paths
	int N{ 1000 };

	// loop over calculations
	for (int i{ 0 }; i < M; i++) {
	
		// store result from each  monte carlo calculation
		std::vector<double> samples;
	
		// loop over different K
		for (int K{ 1 }; K <= 100; K += 10) {

			// calculate and store option value 
			samples.push_back(value_Asian_call(initial_share_price, interest_rate, dividend_rate, volatility, expiration, N, K));

			// store N and K
			if (i == 0) {
				parameter_store.push_back(K);
				value_store.push_back(samples[i]);
			}
		}

		// calculate the mean
		double sum_mean{ 0 };
		for (int i{ 0 }; i < samples.size(); i++) sum_mean += samples[i];
		double sample_mean = sum_mean / M;

		// calculate the variance
		double sum_var{ 0 };
		for (int i{ 0 }; i < samples.size(); i++) sum_var += pow(samples[i] - sample_mean, 2);
		double sample_variance = sum_var / (M - 1.);

		// calculate population mean and std
		double pop_mean = sample_mean;
		double pop_std = sqrt(sample_variance / M);

		// calculate confidence intervals
		upper_confidence.push_back(pop_mean + 2 * pop_std);
		lower_confidence.push_back(pop_mean - 2 * pop_std);
	}

	// open a file stream for writing
	std::ofstream output;

	// open the csv file
	output.open("path dependence with confidence interval K.csv");

	// if the file is open
	if (output.is_open()) {

		// loop over data containers
		for (int i{ 0 }; i < parameter_store.size(); i++) {

			// write data to file
			output << parameter_store[i] << "," << value_store[i] << "," << upper_confidence[i] << "," << lower_confidence[i] << std::endl;
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