#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <string> 
#include <complex>
//#include <stdlib.h>
#include <sstream>

using namespace std;

vector<vector<double>> polar_statevector(vector<std::complex<double>> cartesian) {
	// polar formula for complex numbers from cartesian
	// it translates an entire statevector, so 2^N entries
	// each of them will be a pair (magnitude, argument)
	// declare variables
	float tol = .01;
	int n = cartesian.size();
	double mag, arg, re, im; // auxiliary for each entry
	vector<vector<double>> polar_statevector;

	// loop through state vector
	for (int i = 0; i < n; i++) {
		// get faulty results
		re = real(cartesian[i]);
		im = imag(cartesian[i]);
		mag = sqrt(re * re + im * im);
		if (abs(re) < tol) {
			if (abs(im) < tol) {
				arg = 0;
			}
			else if (im > tol) {
				arg = M_PI / 2;
			}
			else {
				arg = -M_PI / 2;
			}
		}
		else if (im > tol) {
			arg = atan(im / re);
		}
		else {
			arg = M_PI-atan(im / re);
		}

		// recover the ideal solutions
		if (abs(mag - 1) < tol) {
			mag = 1;
		}
		if (abs(sqrt(2) * mag - 1) < tol) {
			mag = 1 / sqrt(2);
		}
		if (abs(2 * mag - 1) < tol) {
			mag = 0.5;
		}
		if (abs(2 * sqrt(2) * mag - 1) < tol) {
			mag = 0.5 / sqrt(2);
		}

		if (abs(arg) < tol) {
			arg = 0;
		}
		if (abs(arg - M_PI) < tol) {
			arg = M_PI;
		}
		vector<double> aux;
		aux.push_back(mag);
		aux.push_back(arg);

		// add entry to polar state vector
		polar_statevector.push_back(aux);
	}

	return polar_statevector;
}

vector<vector<vector<double>>> start_level(int level) {
	// gives initial state and goal state for each level
	// returns states, with states[0] initial state and state[1] goal state

	// initialize variables
	vector<vector<vector<double>>> states;
	vector<vector<double>> initial, goal;
	double aux = 1 / sqrt(2);

	switch (level)	// Initial and goal states are hard coded for each level
	{				// This ain't pretty
	case 0:
		for (int i = 0; i < 8; i++) {
			vector<double> polar1, polar2;
			if (i == 0) {
				polar1.push_back(1.0);
			}
			else {
				polar1.push_back(0.0);
			}
			polar1.push_back(0.0);

			if (i == 6) {
				polar2.push_back(1.0);
				polar2.push_back(M_PI);
			}
			else {
				polar2.push_back(0.0);
				polar2.push_back(0.0);
			}
			initial.push_back(polar1);
			goal.push_back(polar2);
		}
		break;

	case 1:
		for (int i = 0; i < 8; i++) {
			vector<double> polar1, polar2;
			polar1.push_back(0.5 * aux);
			polar1.push_back(0.0);

			if (i == 6) {
				polar2.push_back(aux);
				polar2.push_back(0.0);
			}
			else if (i == 7) {
				polar2.push_back(aux);
				polar2.push_back(M_PI);
			}
			else {
				polar2.push_back(0.0);
				polar2.push_back(0.0);
			}
			initial.push_back(polar1);
			goal.push_back(polar2);
		}
		break;

	case 2:
		for (int i = 0; i < 8; i++) {
			vector<double> polar1, polar2;
			if (i == 0 || i == 3) {
				polar1.push_back(aux);
			}
			else {
				polar1.push_back(0.0);
			}
			polar1.push_back(0.0);

			if (i == 3) {
				polar2.push_back(aux);
				polar2.push_back(-0.5 * M_PI);
			}
			else if (i == 5) {
				polar2.push_back(aux);
				polar2.push_back(0.5 * M_PI);
			}
			else {
				polar2.push_back(0.0);
				polar2.push_back(0.0);
			}
			initial.push_back(polar1);
			goal.push_back(polar2);
		}
		break;
	}

	states.push_back(initial);
	states.push_back(goal);

	return states;
}

vector<string> get_secret_gates(int level) {
	// let the game manager know how to generate the initial state
	// again, this is hard coded so not very beautiful
	
	vector<string> gates;

	switch (level) {
	case 1:
		gates.push_back("h0");
		gates.push_back("h1");
		gates.push_back("h2");
		break;

	case 2:
		gates.push_back("h1");
		gates.push_back("c12");
		break;

	}

	return gates;
}

vector<std::complex<double>> get_end_state(int level) {
	// let the game manager know what the goal state is
	// the end state for each level is hard coded
	vector<std::complex<double>> end_state;
	double aux = 1 / sqrt(2);

	switch (level) {
	case 0:
		for (int i = 0; i < 8; i++) {
			if (i == 6) {
				std::complex<double> e(-1.0, 0.0);
				end_state.push_back(e);
			}
			else {
				std::complex<double> e(0.0, 0.0);
				end_state.push_back(e);
			}
		}
		break;

	case 1:
		for (int i = 0; i < 8; i++) {
			if (i == 6) {
				std::complex<double> e(aux, 0.0);
				end_state.push_back(e);
			}
			else if (i == 7) {
				std::complex<double> e(-aux, 0.0);
				end_state.push_back(e);
			}
			else {
				std::complex<double> e(0.0, 0.0);
				end_state.push_back(e);
			}
		}
		break;

	case 2:
		for (int i = 0; i < 8; i++) {
			if (i == 3) {
				std::complex<double> e(0.0, -aux);
				end_state.push_back(e);
			}
			else if (i == 5) {
				std::complex<double> e(0.0, aux);
				end_state.push_back(e);
			}
			else {
				std::complex<double> e(0.0, 0.0);
				end_state.push_back(e);
			}
		}
		break;
	}

	return end_state;
}

int mypow(int x, int p)
{
	if (p == 0) return 1;
	if (p == 1) return x;

	int tmp = mypow(x, p / 2);
	if (p % 2 == 0) return tmp * tmp;
	else return x * tmp * tmp;
}

vector<int> x_permut(int n, int m) {
	// performs permutation equivalent to X gate
	// args:
	//	n (int): number of qubits
	//	m (int): qubit on which gate is applied
	// rets:
	//	perm (vector int): permutation of 2**n elements

	// initialize variables
	vector<int> perm;
	for (int i = 0; i < mypow(2, n); i++) {
		perm.push_back(i);
	}

	// loop through entries
	// avoid matrix multiplication, quite fancy footwork
	for (int i = 0; i < mypow(2, n - 1); i++) {
		int I = 2 * i - i % mypow(2,m);
		int J = I + mypow(2, m);
		int a = perm[I];
		perm[I] = perm[J];
		perm[J] = a;
	}

	return perm;
}

vector<int> cx_permut(int n, int c, int t) {
	// performs permutation equivalent to CX gate
	// args:
	//	n (int): number of qubits
	//	c (int): control qubit
	//	t (int): target qubit
	// rets:
	//	perm (vector int): permutation of 2**n elements
	vector<int> perm;
	for (int i = 0; i < mypow(2, n); i++) {
		perm.push_back(i);
	}

	// cf. x_permut
	for (int i = 0; i < mypow(2, n - 2); i++) {
		int I = mypow(2, c) + i % mypow(2, c) + ((i - i % mypow(2, c)) * 2) % mypow(2, t) + 2 * ((i - i % mypow(2, c)) * 2 - ((2 * (i - i % mypow(2, c)) % mypow(2, t))));
		int J = I + mypow(2, t);
		int a = perm[I];
		perm[I] = perm[J];
		perm[J] = a;
	}

	return perm;
}