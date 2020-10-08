#include <iostream>
#include "microqiskit.h"
#include <math.h>
#include <vector>
//#include <string> 
//#include <stdlib.h>
//#include <sstream>

using namespace std;

class Interpreter {

public:

	int n_qubits;
	vector<string> gate_list;
	QuantumCircuit qc;

	void set_qubits(int n) {
		// setup qubits
		n_qubits = n;
	}

	Simulator make_circuit(vector<string> gate_list) {
		// make a quantum circuit
		qc.set_registers(n_qubits);

		// initialize variables
		vector<double> real_ket;
		// 2 ^ n_qubits
		int counts = pow(2, n_qubits);

		// build initial vector
		for (int i = 0; i < counts; i = i + 1) {
			if (i == 0) {
				real_ket.push_back(1.0);
			}
			else {
				real_ket.push_back(0.0);
			}
		}
		qc.initialize(real_ket);

		// add gates to quantum circuit
		for (int i = 0; i < gate_list.size(); i = i + 1) {
			string gate_str = gate_list[i];

			if (gate_str[0] == 'x') {
				qc.x(int(gate_str[1]) - 48);
			}
			else if (gate_str[0] == 'y') {
				qc.y(int(gate_str[1]) - 48);
			}
			else if (gate_str[0] == 'z') {
				qc.z(int(gate_str[1]) - 48);
			}
			else if (gate_str[0] == 'h') {
				qc.h(int(gate_str[1]) - 48);
			}
			else if (gate_str[0] == 's') {
				qc.rz(M_PI / 2, int(gate_str[1]) - 48);
			}
			else if (gate_str[0] == 't') {
				qc.rz(M_PI / 4, int(gate_str[1]) - 48);
			}
			else if (gate_str[0] == 'c') {
				qc.cx(int(gate_str[1]) - 48, int(gate_str[2]) - 48);
			}
		}

		Simulator result(qc);
		return result;

	}

	vector<string> get_statevector_string() {
		// get statevectors 
		Simulator result(qc);

		vector<std::complex<double>> ket = result.get_statevector();

		// write vectors to string array
		vector<string> vectors;
		for (int j = 0; j < ket.size(); j++) {
			vectors.push_back("(" + std::to_string(real(ket[j])) + ") + (" + std::to_string(imag(ket[j])) + ")" + "*i");
		}

		return vectors;
	}


};

vector<std::complex<double>> quick_states(int n_qubits, vector<string> gate_list) {
	// initialize class
	Interpreter intp;
	// set number of qubits
	intp.set_qubits(n_qubits);
	// setup circuit
	Simulator result = intp.make_circuit(gate_list);
	return result.get_statevector();
}