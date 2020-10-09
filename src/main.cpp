#include "interpreter.h"
#include <vector>
#include <iostream>


int main() {

	cout << "\n= 1 =" << endl;

	// initialize class
	Interpreter intp;
	// set number of qubits
	intp.set_qubits(2);

	// make gates string
	vector<string> gates;
	gates.push_back("h0");
	//gates.push_back("c01");

	cout << "\n= 2 =" << endl;

	// setup circuit
	Simulator results = intp.make_circuit(gates);

	cout << "\n= 3 =" << endl;

	vector<std::complex<double>> vectors = results.get_statevector();

	cout << "\n= 4 =" << endl;

	cout << "Statevectors:" << endl;
	for (int j = 0; j < vectors.size(); j++) {
		cout << vectors[j] << endl;
	}

	return 0;
} // clang++-7 -pthread -std=c++17 -o main main2.cpp
