#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <string> 
#include <complex>
//#include <stdlib.h>
#include <sstream>

using namespace std;

vector<double> polar_statevector(vector<std::complex<double>> cartesian) {
	// polar formula for complex numbers from cartesian

	// declare variables
	float tol = .01;
	int n = cartesian.size();
	double mag, arg, re, im; // auxiliary for each entry
	vector<double> polar_statevector;

	// loop through state vector
	for (int i = 0; i < n; i++) {
		// get faulty results
		re = real(cartesian[i]);
		im = imag(cartesian[i]);
		mag = re * re + im * im;
		if (abs(re) < tol) {
			if (im > 0) {
				arg = M_PI / 2;
			}
			else if (im < 0) {
				arg = -M_PI / 2;
			}
		}
		else {
			arg = atan(im / re);
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
		vector<double> aux{ mag, arg };

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

	case 1:
		double aux = 1 / sqrt(2);
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

	case 2:
		double aux = 1 / sqrt(2);
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

	case 2:
		gates.push_back("h1");
		gates.push_back("c12");

	}

	return gates;
}

vector<std::complex<double>> get_end_state(int level) {
	// let the game manager know what the goal state is
	// the end state for each level is hard coded
	vector<std::complex<double>> end_state;

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

	case 1:
		double aux = 1 / sqrt(2);
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

	case 2:
		double aux = 1 / sqrt(2);
		for (int i = 0; i < 8, i++) {
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
	}

	return end_state;
}