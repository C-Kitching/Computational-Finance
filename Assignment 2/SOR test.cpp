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
double crank_nicolson(const double& T, const double& r, const double& S0, const double& X, const double& sigma, const int& i_max, const int& j_max);


// Begin main program
int main()
{
	// declare and initialise parameters
	double T{ 1 };
	double r{ 0.05 };
	double S0{ 1 };
	double X{ 2 };
	double sigma{ 0.4 };

	// declare and initialise grid parameters
	int i_max{ 4 };
	int j_max{ 4 };

	auto start = std::chrono::steady_clock::now();  // get start time

	// get option value
	double option_value = crank_nicolson(T, r, S0, X, sigma, i_max, j_max);

	auto finish = std::chrono::steady_clock::now();  // get finish time
	auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>> (finish - start);  // time difference

	// output calculation time
	std::cout << elapsed.count() << std::endl;

	// output result
	std::cout << "V(S = " << S0 << ", t = 0) = " << std::setprecision(10) << option_value << std::endl;



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

	if (iterations == max_iter) std::cout << "\n NOT CONVERGED \n";
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
double crank_nicolson(const double& T,  const double& r, const double& S0, const double& X, const double& sigma, const int& i_max, const int& j_max)
{
	// declare and initialise local parameters (dS, dt)
	double S_max = 2 * X;
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
	double tolerence = 1e-8;
	double omega = 1;

	// initialise our stock prices
	for (int j{ 0 }; j <= j_max; j++) S[j] = j * dS;

	// initialise final conditions on the option price
	for (int j{ 0 }; j <= j_max; j++) {
		v_old[j] = std::max(X - S[j], 0.);
		v_new[j] = std::max(X - S[j], 0.);
	}

	// loop over the time levels
	for (int i{ i_max - 1 }; i >= 0; i--) {

		// declare matrix parameter vectors
		std::vector<double> a(j_max + 1);
		std::vector<double> b(j_max + 1);
		std::vector<double> c(j_max + 1);
		std::vector<double> d(j_max + 1);

		// initial values at j = 0
		a[0] = 0.;
		b[0] = 1.;
		c[0] = 0.;
		d[0] = X * exp(-r * (T - i * dt));

		// loop through middling values of j
		for (int j{ 1 }; j <= j_max - 1; j++) {

			a[j] = 0.25 * (sigma * sigma * j * j - r * j);
			b[j] = -0.5 * sigma * sigma * j * j - 0.5 * r - 1. / dt;
			c[j] = 0.25 * (sigma * sigma * j * j + r * j);
			d[j] = -a[j] * v_old[j - 1] - (b[j] + 2. / dt) * v_old[j] - c[j] * v_old[j + 1];
		}

		// initialise values at j = j_max
		a[j_max] = 0.;
		b[j_max] = 1.;
		c[j_max] = 0.;
		d[j_max] = 0.;

		for (int j = 0; j <= j_max; j++)
		{
			std::cout << j << " " << a[j] << " " << b[j] << " " << c[j] << " " << d[j] << std::endl;
		}

		// solve
		iterations = 0;
		SOR_solve(a, b, c, d, v_new, max_iter, tolerence, omega, iterations);
		std::cout << " Solved after " << iterations << " iterations.\n";

		std::cout << "vold" << std::endl;
		for (int i{ 0 }; i < v_old.size(); i++) {
			std::cout << i << "," << v_old[i] << std::endl;
		}

		// set old to new
		v_old = v_new;

		std::cout << "vnew" << std::endl;
		for (int i{ 0 }; i < v_old.size(); i++) {
			std::cout << i << "," << v_old[i] << std::endl;
		}
	}

	double optionValue;
	{
		int jStar = S0 / dS;
		double sum = 0.;
		sum += (S0 - S[jStar]) / dS * v_new[jStar + 1];
		sum += (S[jStar + 1] - S0) / dS * v_new[jStar];
		optionValue = sum;
	}

	// use lagrange interpolation to get estimated option value
	double optionValue2 = lagrange_interpolation(v_new, S, S0, 16);

	return optionValue2;
}