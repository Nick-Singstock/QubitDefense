#include <iostream>
#include "microqiskit.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <string> 
//#include <stdlib.h>
#include <sstream>

using namespace std;

vector<double> polar_statevector(vector<std::complex<double>> cartesian) {
	// polar formula for complex numbers from cartesian

	// declare variables
	float tol = .1;
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
		vector<double> aux {mag, arg}

		// add entry to polar state vector
		polar_statevector.push_back(aux)
	}

	return polar_statevector;
}

