#include <iostream>
#include "microqiskit.h"
#include <math.h>
#include <vector>
#include "game_state.h"

using namespace std;

// make variables global
//int level;
//vector<string> gate_list;

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
		for (int i = 0; i < gate_list.size(); i++) {
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

bool is_equal(double s1, double s2) {
	double threshold = 0.01;
	if (s1 - s2 < threshold && s1 - s2 > -threshold) {
		return true;
	}
	return false;
}

bool check_states(vector<std::complex<double>> current_state, vector<std::complex<double>> end_state) {
	for (int i = 0; i < current_state.size(); i++) {
		if (is_equal(real(current_state[i]), real(end_state[i])) == false) {
			return false;
		}
		else if (is_equal(imag(current_state[i]), imag(end_state[i])) == false) {
			return false;
		}
	}
	return true;
}

// function to combine list of logic gate operation strings
vector<string> combine_gates(vector<string> list1, vector<string> list2) {
	vector<string> new_list;
	for (int i = 0; i < list1.size(); i++) {
		new_list.push_back(list1[i]);
	}
	for (int i = 0; i < list2.size(); i++) {
		new_list.push_back(list2[i]);
	}
	return new_list;
}


// will manage the current status of the game
class game_manager {
public:
	int n_qubits = 3; // hard coded for now but could come from level

	int level = 1;
	vector<string> gate_list;

	// get secret states and end state goal, as defined by the level
	vector<string> secret_gates = get_secret_gates(level);
	vector<std::complex<double>> end_state = get_end_state(level);

	// combine secret_gates and gate_list to get full_gate_list
	vector<string> full_gate_list = combine_gates(secret_gates, gate_list); 

	// get current statevector
	vector<std::complex<double>> current_state = quick_states(n_qubits, full_gate_list);
	
	// call to determine if game is won
	bool game_won() {
		// get current statevector
		vector<std::complex<double>> current_state = quick_states(n_qubits, full_gate_list);
		return check_states(current_state, end_state);
	}

	void update_level(int new_level) {
		level = new_level;
		// get secret states and end state goal, as defined by the level
		vector<string> secret_gates = get_secret_gates(level);
		vector<std::complex<double>> end_state = get_end_state(level);
	}

	void add_gate(string new_gate) {
		gate_list.push_back(new_gate);
		// combine secret_gates and gate_list to get full_gate_list
		vector<string> full_gate_list = combine_gates(secret_gates, gate_list);
	}

	void update_gates(vector<string> new_gate_list) {
		gate_list = new_gate_list;
		// combine secret_gates and gate_list to get full_gate_list
		vector<string> full_gate_list = combine_gates(secret_gates, gate_list);
	}

	// call to get current statevector of game
	vector<vector<double>> statevector() {
		// get current statevector
		vector<std::complex<double>> current_state = quick_states(n_qubits, full_gate_list);
		return polar_statevector(current_state);
	}

	// call to get goal state
	vector<vector<double>> goal_state() {
		return polar_statevector(end_state);
	}
	// call to get complex goal state for testing purposes
	vector<std::complex<double>> goal_state_complex() {
		return end_state;
	}

};
