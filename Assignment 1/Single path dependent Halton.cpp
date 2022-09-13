// HEADER
// Student ID: 10134521
// Title: Assignment 1
// Date Created: 18/03/21
// Last Edited:


// math constants
#define _USE_MATH_DEFINES


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

// generate first N prime numbers
std::vector<int> prime(const int& N);

// generate Halton sequence
std::vector<double> Halton_sequence(const int& basis, const int& size);

// value Asian call
double value_Asian_call(const double& initial_share_price, const double& interest_rate, const double& dividend_rate, const double& volatility,
	const double& expiration, const int& N, const int& K);


// Begin main program
int main()
{
	// define parameters
	double expiration{ 1.25 };
	double volatility{ 0.37 };
	double interest_rate{ 0.03 };
	double dividend_rate{ 0.04 };
	double initial_share_price{ 900 };
	int K{ 35 };  // points in the sample path
	int N{ 2000 };  // number of MC paths

	auto start = std::chrono::steady_clock::now();  // get start time

	// value the option
	double value = value_Asian_call(initial_share_price, interest_rate, dividend_rate, volatility, expiration, N, K);

	auto finish = std::chrono::steady_clock::now();  // get finish time

	auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>> (finish - start);  // convert into seconds

	std::cout << "Time = " << elapsed.count() << std::endl;  // output time

	// output results
	std::cout << value << std::endl;

	return 0;
}  // End main progrma


// Function definitions

// generate first N prime numbers
std::vector<int> prime(const int& N) 
{
	// initialise
	int test_number{ 2 };
	std::vector<int> primes;
	
	// while we havent generated N primes
	while (primes.size() < N) {
	
		// assume number is prime
		int flag{ 1 };

		// for loop over all possible factors
		for (int j{ 2 }; j <= test_number / 2; j++) {
		
			// if we find a factor
			if (test_number % j == 0) {
				flag = 0;  // not a prime
				break;
			}
		}
	
		// if number is prime, store it
		if (flag == 1) primes.push_back(test_number);
	
		// update test number
		test_number++;
	}

	return primes;
}

// generate Halton sequence
std::vector<double> Halton_sequence(const int& basis, const int& size)
{
	// declare vector to return
	std::vector<double> Halton;

	// generate vector of size N
	for (int i{ 1 }; i <= size; i++) {

		// initialise variables
		double temp{ 1 };
		double Halton_number{ 0 };
		int index{ i };

		// calculate Halton number at index
		while (index > 0) {

			temp /= basis;
			Halton_number += temp * (index % basis);
			index /= basis;
		}

		// record the number
		Halton.push_back(Halton_number);
	}

	return Halton;
}

// value Asian call
double value_Asian_call(const double& initial_share_price, const double& interest_rate, const double& dividend_rate, const double& volatility,
	const double& expiration, const int& N, const int& K)
{
	// declare random number generator
	static std::mt19937 rnd;

	// declare the normal distrubtion
	std::normal_distribution<double> ND(0., 1.);

	// set the basis
	int basis_1{ 2 };
	int basis_2{ 3 };

	// generate Halton sequence of length N*K
	std::vector<double> Halton_1 = Halton_sequence(basis_1, N * K);
	std::vector<double> Halton_2 = Halton_sequence(basis_2, N * K);

	// convert to random normal with Box-Muller
	std::vector<double> random_1;
	std::vector<double> random_2;
	for (int i{ 0 }; i < N * K; i++) {
		random_1.push_back(cos(2 * M_PI * Halton_2[i]) * pow(-2 * log(Halton_1[i]), 0.5));
		random_2.push_back(sin(2 * M_PI * Halton_1[i]) * pow(-2 * log(Halton_2[i]), 0.5));
	}

	// initalise sum to zero
	double sum{ 0 };

	// counter to use correct random number
	int counter_1{ 0 };
	int counter_2{ 0 };

	//pseudorandom number stores
	std::vector<double> pseudo_1;
	std::vector<double> pseudo_2;

	// halton random stores
	std::vector<double> halton_1;
	std::vector<double> halton_2;

	// final path values 
	std::vector<double> ST1;
	std::vector<double> ST2;
	std::vector<double> ST3;
	std::vector<double> ST4;

	// loop over all MC paths
	for (int i{ 0 }; i < N; i++) {

		// time interval
		double dt{ expiration / K };

		// initialise stock paths
		std::vector<double> stock_path_1;
		std::vector<double> stock_path_2;
		std::vector<double> stock_path_3;
		std::vector<double> stock_path_4;

		// append initial value
		stock_path_1.push_back(initial_share_price);
		stock_path_2.push_back(initial_share_price);
		stock_path_3.push_back(initial_share_price);
		stock_path_4.push_back(initial_share_price);

		// generate paths
		for (int j{ 1 }; j <= K; j++) {

			// generate random number
			double phi1 = random_1[counter_1];
			double phi2 = random_2[counter_2];
			halton_1.push_back(phi1);
			halton_2.push_back(phi2);

			// generate and store pseudorandom numbers 
			double phi3 = ND(rnd);
			double phi4 = ND(rnd);
			pseudo_1.push_back(phi3);
			pseudo_2.push_back(phi4);

			// gemerate stock path with Halton
			stock_path_1.push_back(stock_path_1[j - 1] * exp((interest_rate - dividend_rate -
				0.5 * pow(volatility, 2)) * dt + volatility * phi1 * pow(dt, 0.5)));
			stock_path_2.push_back(stock_path_2[j - 1] * exp((interest_rate - dividend_rate -
				0.5 * pow(volatility, 2)) * dt + volatility * phi2 * pow(dt, 0.5)));

			// gemerate stock path with pseduo
			stock_path_3.push_back(stock_path_1[j - 1] * exp((interest_rate - dividend_rate -
				0.5 * pow(volatility, 2)) * dt + volatility * phi3 * pow(dt, 0.5)));
			stock_path_4.push_back(stock_path_2[j - 1] * exp((interest_rate - dividend_rate -
				0.5 * pow(volatility, 2)) * dt + volatility * phi4 * pow(dt, 0.5)));

			// increment counters
			counter_1++;
			counter_2++;
		}

		// calculate A
		double A1, A2;
		double A1_sum{ 0 }, A2_sum{ 0 };
		for (int i{ 1 }; i <= K; i++) {
			A1_sum += stock_path_1[i];
			A2_sum += stock_path_2[i];
		}
		A1 = A1_sum / K;
		A2 = A2_sum / K;

		// add in the payoff
		sum += std::max(stock_path_1.back() - A1, 0.);
		sum += std::max(stock_path_2.back() - A2, 0.);

		ST1.push_back(stock_path_1.back());
		ST2.push_back(stock_path_2.back());
		ST3.push_back(stock_path_3.back());
		ST4.push_back(stock_path_4.back());

	}

	// open a file stream for writing
	std::ofstream output;

	// open the csv file
	output.open("normal_test.csv");

	// if the file is open
	if (output.is_open()) {

		for (int i{ 0 }; i < random_1.size(); i++) {
			output << pseudo_1[i] << "," << pseudo_2[i] << "," << halton_1[i] << "," << halton_2[i] 
				<< "," << ST1[i] << "," << ST2[i] << "," << ST3[i] << "," << ST4[i] << std::endl;
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

	// average over all paths
	return exp(-interest_rate * expiration) * sum / (2.*N);
}