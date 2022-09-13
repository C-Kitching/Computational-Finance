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


// Function declerations

// value Asian call option
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

	int M{ 100 };  // number of calculations
	std::vector<double> samples;  // store result from each calculation

	int N{ 788096 };  // number of monte carlo simulations to perform

	// loop over the number of calculations
	for (int i{ 0 }; i < M; i++) {

		// store monte carlo result
		samples.push_back(value_Asian_call(initial_share_price, interest_rate, dividend_rate, volatility, expiration, N, K));
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
	double upper_95 = pop_mean + 2 * pop_std;
	double lower_95 = pop_mean - 2 * pop_std;

	// output results
	std::cout << "95% confidence result is in [" << lower_95 << "," << upper_95 << "] with N = " << N << std::endl;
	std::cout << "V = " << pop_mean << std::endl;

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
	return exp(-interest_rate * expiration) * sum / (2. * N);
}