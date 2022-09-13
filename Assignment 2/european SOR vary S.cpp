// Includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <chrono>


// Function declerations

// calculate theta(t)
double theta(const double& mu, const double& X, const double& dt, const int& i);

// SOR solver
void SOR_solve(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d,
	std::vector<double>& solution, const int& max_iter, const double& tolerence, const double& omega, int& iterations);

// generic lagrange interpolation
double lagrange_interpolation(const std::vector<double>& y, const std::vector<double>& x, double x0, unsigned int n);

// Crank Nicolson Finite Difference
double crank_nicolson(const double& T, const double& F, const double& R, const double& r, const double& kappa, const double& mu,
	const double& S0, const double& X, const double& C, const double& alpha, const double& beta, const double& sigma, const int& i_max, const int& j_max);

// calculate d1
double d1(const double& share_price, const double& strike_price, const double& interest_rate,
	const double& volatility, const double& expiration, const double& time);

// calculate d2
double d2(const double& share_price, const double& strike_price, const double& interest_rate,
	const double& volatility, const double& expiration, const double& time);

// analystical call
double analytic_call(const double& share_price, const double& strike_price, const double& interest_rate,
	const double& volatility, const double& expiration, const double& time);

// normal cummulative distribution
double norm_cumm(const double& x);


// Begin main program
int main()
{
	// declare and initialise parameters
	double T{ 2 };
	double F{ 50 };
	double R{ 1 };
	double r{ 0.0114 };
	double kappa{ 0.125 };
	double mu{ 0.0174 };
	double S0{ 40 };
	double X{ 50.5 };
	double C{ 0.285 };
	double alpha{ 0.01 };
	double beta{ 0.869 };
	double sigma{ 0.668 };

	// test variables
	kappa = 0;
	beta = 1;
	C = 0;

	// maximum S
	double S_max = 4 * X;

	// declare and initialise grid parameters
	int i_max{ 10 };
	int j_max{ 10 };

	// containers for data
	std::vector<double> S_store;
	std::vector<double> V_store;
	std::vector<double> analytic_store;

	// number of data points 
	int n{ 100 };

	auto start = std::chrono::steady_clock::now();  // get start time

	// loop over data points
	for (int S{ 0 }; S <= S_max; S += S_max / n) {

		S_store.push_back(S);
		V_store.push_back(crank_nicolson(T, F, R, r, kappa, mu, S, X, C, alpha, beta, sigma, i_max, j_max));
		analytic_store.push_back(F*exp(-r*T) + analytic_call(S, F, r, sigma, T, 0));
	}

	auto finish = std::chrono::steady_clock::now();  // get finish time
	auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>> (finish - start);  // time difference

	// output calculation time
	std::cout << elapsed.count() << std::endl;

	// open a file stream for writing
	std::ofstream output;

	// open the csv file
	output.open("european SOR vary S.csv");

	// if the file is open
	if (output.is_open()) {

		// write data to file
		for (int i{ 0 }; i < S_store.size(); i++) {
			output << S_store[i] << "," << V_store[i] << "," << analytic_store[i] << std::endl;
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

}  // End main program


// Function definitions

// calculate theta(t)
double theta(const double& mu, const double& X, const double& dt, const int& i)
{
	return (1 + mu) * X * exp(mu * i * dt);
}

// SOR solver
void SOR_solve(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d,
	std::vector<double>& solution, const int& max_iter, const double& tolerence, const double& omega, int& iterations)
{
	// get size of vector
	int n = solution.size() - 1;

	// decalre y
	double y;

	// SOR loop
	for (iterations; iterations < max_iter; iterations++) {

		// reset error to 0
		double error = 0;

		// initial value
		y = (d[0] - c[0] * solution[1]) / b[0];
		solution[0] = solution[0] + omega * (y - solution[0]);

		// middelling values
		for (int j = 1; j < n; j++)
		{
			double y = (d[j] - a[j] * solution[j - 1] - c[j] * solution[j + 1]) / b[j];
			solution[j] = solution[j] + omega * (y - solution[j]);
		}

		// final value
		y = (d[n] - a[n] * solution[n - 1]) / b[n];
		solution[n] = solution[n] + omega * (y - solution[n]);

		// calculate residual norm ||r|| as sum of absolute values
		error += std::fabs(d[0] - b[0] * solution[0] - c[0] * solution[1]);
		for (int j = 1; j < n; j++) {
			error += std::fabs(d[j] - a[j] * solution[j - 1] - b[j] * solution[j] - c[j] * solution[j + 1]);
		}
		error += std::fabs(d[n] - a[n] * solution[n - 1] - b[n] * solution[n]);

		// make an exit condition when solution found
		if (error < tolerence) return;
	}

	if (iterations = max_iter) std::cout << "No convergence!" << std::endl;
}

// generic lagrange interpolation
double lagrange_interpolation(const std::vector<double>& y, const std::vector<double>& x, double x0, unsigned int n)
{
	if (x.size() < n) return lagrange_interpolation(y, x, x0, x.size());
	if (n == 0)throw;

	// local parameters
	int nHalf = n / 2;
	double dx = x[1] - x[0];

	// calculate j star
	int jStar;
	if (n % 2 == 0) jStar = int((x0 - x[0]) / dx) - (nHalf - 1);  // even degree
	else jStar = int((x0 - x[0]) / dx + 0.5) - (nHalf);  // odd degree

	jStar = std::max(0, jStar);
	jStar = std::min(int(x.size() - n), jStar);

	if (n == 1)return y[jStar];

	double temp = 0.;
	for (unsigned int i = jStar; i < jStar + n; i++) {

		double  int_temp;
		int_temp = y[i];

		for (unsigned int j = jStar; j < jStar + n; j++) {

			if (j == i) { continue; }
			int_temp *= (x0 - x[j]) / (x[i] - x[j]);
		}
		temp += int_temp;
	}  // end of interpolate

	return temp;
}

// Crank Nicolson Finite Difference
double crank_nicolson(const double& T, const double& F, const double& R, const double& r, const double& kappa, const double& mu,
	const double& S0, const double& X, const double& C, const double& alpha, const double& beta, const double& sigma, const int& i_max, const int& j_max)
{
	// declare and initialise local parameters (dS, dt)
	double S_max = 5 * X;
	double dS = S_max / j_max;
	double dt = T / i_max;

	// declare iteration number
	int iterations{ 0 };

	// create storage for stock price and old and new option price
	std::vector<double> S(j_max + 1);
	std::vector<double> v_old(j_max + 1);
	std::vector<double> v_new(j_max + 1);

	// params for SOR solve
	int max_iter = 10000;
	double tolerence = 1e-10;
	double omega = 1;

	// initialise our stock prices
	for (int j{ 0 }; j <= j_max; j++) S[j] = j * dS;

	// initialise final conditions on the option price
	for (int j{ 0 }; j <= j_max; j++) {
		v_old[j] = std::max(F, R * S[j]);
		v_new[j] = std::max(F, R * S[j]);
	}

	// loop over the time levels
	for (int i{ i_max - 1 }; i >= 0; i--) {

		// declare matrix parameter vectors
		std::vector<double> a(j_max + 1);
		std::vector<double> b(j_max + 1);
		std::vector<double> c(j_max + 1);
		std::vector<double> d(j_max + 1);

		// initial values at j = 0
		a[0] = 0;
		b[0] = -(1 / dt) - (kappa * theta(mu, X, dt, i) / dS) - (r / 2);
		c[0] = kappa * theta(mu, X, dt, i);
		d[0] = (-(1 / dt) + (r / 2)) * v_old[0] - C * exp(-i * dt);

		// loop through middling values of j
		for (int j{ 1 }; j <= j_max - 1; j++) {

			a[j] = -0.25 * pow(sigma, 2) * pow(j, 2 * beta) * pow(dS, 2 * (beta - 1)) + (kappa / (4 * dS)) * (theta(mu, X, dt, i) - j * dS);
			b[j] = (1 / dt) + 0.5 * pow(sigma, 2) * pow(j, 2 * beta) * pow(dS, 2 * (beta - 1)) + (r / 2);
			c[j] = -0.25 * pow(sigma, 2) * pow(j, 2 * beta) * pow(dS, 2 * (beta - 1)) - (kappa / (4 * dS)) * (theta(mu, X, dt, i) - j * dS);
			d[j] = (0.25 * pow(sigma, 2) * pow(j, 2 * beta) * pow(dS, 2 * (beta - 1)) - (kappa / (4 * dS)) * (theta(mu, X, dt, i) - j * dS)) * v_old[j - 1]
				+ ((1 / dt) - 0.5 * pow(sigma, 2) * pow(j, 2 * beta) * pow(dS, 2 * (beta - 1)) - (r / 2)) * v_old[j]
				+ (0.25 * pow(sigma, 2) * pow(j, 2 * beta) * pow(dS, 2 * (beta - 1)) + (kappa / (4 * dS)) * (theta(mu, X, dt, i) - j * dS)) * v_old[j + 1]
				+ C * exp(-alpha * i * dt);
		}

		// initialise values at j = j_max
		a[j_max] = 0;
		b[j_max] = 1;
		c[j_max] = 0;
		d[j_max] = R * (S[j_max] - X) * exp(-(kappa + r) * (T - i * dt)) + (X * R + (C / (alpha + r)) * exp(-alpha * i * dt)) * exp(-r * (T - i * dt)) + (C / (alpha + r)) * exp(-alpha * i * dt);

		// solve
		iterations = 0;
		SOR_solve(a, b, c, d, v_new, max_iter, tolerence, omega, iterations);

		// set old to new
		v_old = v_new;
	}

	// use lagrange interpolation to get estimated option value
	double option_value = lagrange_interpolation(v_new, S, S0, 16);
	return option_value;
}

// calculate d1
double d1(const double& share_price, const double& strike_price, const double& interest_rate,
	const double& volatility, const double& expiration, const double& time)
{
	return (log(share_price / strike_price) + (interest_rate + pow(volatility, 2) / 2) * (expiration - time)) / (volatility * pow(expiration - time, 0.5));
}

// calculate d2
double d2(const double& share_price, const double& strike_price, const double& interest_rate,
	const double& volatility, const double& expiration, const double& time)
{
	return d1(share_price, strike_price, interest_rate, volatility, expiration, time) - volatility * pow(expiration - time, 0.5);
}

// analystical call
double analytic_call(const double& share_price, const double& strike_price, const double& interest_rate,
	const double& volatility, const double& expiration, const double& time)
{
	double d1_val = d1(share_price, strike_price, interest_rate, volatility, expiration, time);
	double d2_val = d2(share_price, strike_price, interest_rate, volatility, expiration, time);

	return share_price * norm_cumm(d1_val) - strike_price * exp(-interest_rate * (expiration - time)) * norm_cumm(d2_val);
}

// normal cummulative distribution
double norm_cumm(const double& x)
{
	return 0.5 * erfc(-x / pow(2, 0.5));
}