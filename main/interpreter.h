#include <iostream>
#include "microqiskit.h"
#include <math.h>
#include <vector>
#include <string> 
//#include <stdlib.h>
#include <sstream>

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

	void make_circuit(vector<string> gate_list) {
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
			cout << gate_str << endl;

			if (gate_str[0] == 'x') {
				qc.x(int(gate_str[1]));
			}
			else if (gate_str[0] == 'y') {
				qc.y(int(gate_str[1]));
			}
			else if (gate_str[0] == 'z') {
				qc.z(int(gate_str[1]));
			}
			else if (gate_str[0] == 'h') {
				stringstream nss(gate_str[1]);
				//const char n = gate_str[1];
				//cout << n << endl;
				int x;// = atoi (n);
				nss >> x;
				cout << "Value x = " << x << endl;
				//istringstream(gate_str[1]) >> x;

				//int x = stoi(gate_str[1])
				qc.h(x);
			}
			else if (gate_str[0] == 'c') {
				qc.cx(int(gate_str[1]), int(gate_str[2]));
			}
		}

		Simulator result(qc);
		cout << result.get_qiskit() << endl;

	}

	vector<string> get_statevectors() {
		// get statevectors 
		Simulator result(qc);

		cout << "\n= 3.1 =" << endl;

		vector<std::complex<double>> ket = result.get_statevector();

		cout << "\n= 3.2 =" << endl;

		// write vectors to string array
		vector<string> vectors;
		for (int j = 0; j < ket.size(); j++) {
			vectors.push_back("(" + std::to_string(real(ket[j])) + ") + (" + std::to_string(imag(ket[j])) + ")" + "*i");
		}

		return vectors;
	}


};
