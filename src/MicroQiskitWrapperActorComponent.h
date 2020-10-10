// Fill out your copyright notice in the Description page of Project Settings.

#pragma once
#include <iostream>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <complex>  
#include <iostream>
#include <bitset>
#include <ctime>
#include <map>

#include "CoreMinimal.h"
#include "Components/ActorComponent.h"
#include "Math/IntPoint.h"
#include "MicroQiskitWrapperActorComponent.generated.h"



//#define max(x,y) ((x)>(y)?(x):(y))

using namespace std;


const int n_qubits_max = 20; // limit on qubit number
// TO DO: Remove this! It is only used in conversions to bit strings, because it doesn't like variables.

class QuantumCircuit {

public:

	int n_qubits, n_bits;
	vector<vector<string>> data;

	// TO DO: Could this be done by a constructor, and still be consistent with the usage of QuantumCircuit objects in Simulator?
	void set_registers(int n, int m = 0) {
		n_qubits = n;
		n_bits = m;
		// TO DO: Only the cases n_qubits=n_bits and n_bits=0 are allowed in MicroQiskit.
		// Abort and explain if the user provides other inputs.
	}

	void add(QuantumCircuit qc2) {
		//define max(x, y) ((x) > (y) ? (x) : (y))
		//n_bits = max(n_bits, qc2.n_bits);
		n_bits = ((n_bits) > (qc2.n_bits) ? (n_bits) : (qc2.n_bits));
		for (int g = 0; g < qc2.data.size(); g++) {
			data.push_back(qc2.data[g]);
		}
		// TO DO: It is only possible to add circuits with equal n_qubits in MicroQiskit, and qc2.n_bits cannot be non-zero if qc.n_bits is.
		// Abort and explain if the user provides other inputs.
	}

	// TO DO: Generalize initialize function to complex numbers
	void initialize(vector<double> ket) {
		vector<string> gate;
		gate.push_back("initialize");
		for (int j = 0; j < ket.size(); j++) {
			gate.push_back(to_string(ket[j]));
		}
		data.push_back(gate);
	}
	void x(int q) {
		vector<string> gate;
		gate.push_back("x");
		gate.push_back(to_string(q));
		data.push_back(gate);
	}
	void rx(double theta, int q) {
		vector<string> gate;
		gate.push_back("rx");
		gate.push_back(to_string(theta));
		gate.push_back(to_string(q));
		data.push_back(gate);
	}
	void h(int q) {
		vector<string> gate;
		gate.push_back("h");
		gate.push_back(to_string(q));
		data.push_back(gate);
	}
	void cx(int s, int t) {
		vector<string> gate;
		gate.push_back("cx");
		gate.push_back(to_string(s));
		gate.push_back(to_string(t));
		data.push_back(gate);
	}
	void measure(int q, int b) {
		vector<string> gate;
		gate.push_back("m");
		gate.push_back(to_string(b));
		gate.push_back(to_string(q));
		data.push_back(gate);
		// TO DO: It is only possible to add measure gates of the form measure(j,j) in MicroQiskit.
		// Abort and explain if the user provides other inputs.
	}
	void rz(double theta, int q) {
		h(q);
		rx(theta, q);
		h(q);
	}
	void ry(double theta, int q) {
		rx(M_PI / 2, q);
		h(q);
		rx(theta, q);
		h(q);
		rx(-M_PI / 2, q);
	}
	void z(int q) {
		rz(M_PI, q);
	}
	void y(int q) {
		z(q);
		x(q);
	}

};

class Simulator {
	// Contains methods required to simulate a circuit and provide the desired outputs.

	vector<vector<double>> simulate(/*QuantumCircuit qc*/) {

		vector<vector<double>> ket;

		for (int j = 0; j < pow(2, qc.n_qubits); j++) {
			vector<double> e;
			for (int k = 0; k <= 2; k++) {
				e.push_back(0.0);
			}
			ket.push_back(e);
		}
		ket[0][0] = 1.0;


		for (int g = 0; g < qc.data.size(); g++) {

			if (qc.data[g][0] == "initialize") {

				for (int j = 0; j < pow(2, qc.n_qubits); j++) {
					ket[j][0] = stof(qc.data[g][j + 1]);
				}

			}
			else if ((qc.data[g][0] == "x") || (qc.data[g][0] == "rx") || (qc.data[g][0] == "h")) {

				int q;
				q = stoi(qc.data[g][qc.data[g].size() - 1]);

				for (int i0 = 0; i0 < pow(2, q); i0++) {
					for (int i1 = 0; i1 < pow(2, qc.n_qubits - q - 1); i1++) {
						int b0, b1;
						b0 = i0 + int(pow(2, q + 1)) * i1;
						b1 = b0 + int(pow(2, q));

						vector<double> e0, e1;
						e0 = ket[b0];
						e1 = ket[b1];

						if (qc.data[g][0] == "x") {
							ket[b0] = e1;
							ket[b1] = e0;
						}
						else if (qc.data[g][0] == "rx") {
							double theta = stof(qc.data[g][1]);
							ket[b0][0] = e0[0] * cos(theta / 2) + e1[1] * sin(theta / 2);
							ket[b0][1] = e0[1] * cos(theta / 2) - e1[0] * sin(theta / 2);
							ket[b1][0] = e1[0] * cos(theta / 2) + e0[1] * sin(theta / 2);
							ket[b1][1] = e1[1] * cos(theta / 2) - e0[0] * sin(theta / 2);
						}
						else if (qc.data[g][0] == "h") {
							for (int k = 0; k <= 2; k++) {
								ket[b0][k] = (e0[k] + e1[k]) / sqrt(2);
								ket[b1][k] = (e0[k] - e1[k]) / sqrt(2);
							}
						}

					}
				}

			}
			else if (qc.data[g][0] == "cx") {
				int s, t, l, h;
				s = stoi(qc.data[g][qc.data[g].size() - 2]);
				t = stoi(qc.data[g][qc.data[g].size() - 1]);
				if (s > t) {
					h = s;
					l = t;
				}
				else {
					h = t;
					l = s;
				}

				for (int i0 = 0; i0 < pow(2, l); i0++) {
					for (int i1 = 0; i1 < pow(2, h - l - 1); i1++) {
						for (int i2 = 0; i2 < pow(2, qc.n_qubits - h - 1); i2++) {

							int b0, b1;
							b0 = i0 + pow(2, l + 1) * i1 + pow(2, h + 1) * i2 + pow(2, s);
							b1 = b0 + pow(2, t);

							vector<double> e0, e1;
							e0 = ket[b0];
							e1 = ket[b1];

							ket[b0] = e1;
							ket[b1] = e0;

						}
					}
				}

			}

		}

		return ket;

	}

	vector<double> get_probs(/*QuantumCircuit qc*/) {

		// TO DO: For get_counts and get_memory (both of which call this function) the circuit should have a full set of measure gates.
		// Abort and explain if the user does not provide this input.

		vector<vector<double>> ket;
		ket = simulate(/*qc*/);

		vector<double> probs;
		for (int j = 0; j < ket.size(); j++) {

			probs.push_back(pow(ket[j][0], 2) + pow(ket[j][1], 2));

		}

		return probs;

	}

public:

	QuantumCircuit qc;
	int shots;

	Simulator(QuantumCircuit qc_in, int shots_in = 1024) {
		srand((unsigned)time(0));
		qc = qc_in;
		shots = shots_in;
	}

	vector<std::complex<double>> get_statevector() {

		vector<vector<double>> ket;
		ket = simulate(/*qc*/);

		vector<std::complex<double>> complex_ket;
		for (int j = 0; j < ket.size(); j++) {

			std::complex<double> e(ket[j][0], ket[j][1]);
			complex_ket.push_back(e);

		}

		return complex_ket;

	}

	vector<string> get_memory() {

		vector<double> probs;
		probs = get_probs(/*qc*/);

		vector<string> memory;

		for (int s = 0; s < shots; s++) {

			double cumu = 0;
			bool un = true;
			double r = double(rand()) / RAND_MAX;

			for (int j = 0; j < probs.size(); j++) {
				cumu += probs[j];
				if ((r < cumu) && un) {
					std::string long_out = std::bitset<n_qubits_max>(j).to_string();
					std::string out = long_out.substr(n_qubits_max - qc.n_qubits, n_qubits_max);
					memory.push_back(out);
					un = false;
				}
			}

		}

		return memory;

	}

	std::map<std::string, int> get_counts() {

		std::map<std::string, int> counts;

		vector<string> memory = get_memory();

		for (int s = 0; s < shots; s++) {
			counts[memory[s]] += 1;
		}

		return counts;

	}

	string get_qiskit() {

		string qiskit_py;

		if (qc.n_bits == 0) {
			qiskit_py += "qc = QuantumCircuit(" + to_string(qc.n_qubits) + ")\n";
		}
		else {
			qiskit_py += "qc = QuantumCircuit(" + to_string(qc.n_qubits) + "," + to_string(qc.n_bits) + ")\n";
		}

		for (int g = 0; g < qc.data.size(); g++) {
			if (qc.data[g][0] == "x") {
				qiskit_py += "qc.x(" + qc.data[g][1] + ")\n";
			}
			else if (qc.data[g][0] == "rx") {
				qiskit_py += "qc.rx(" + qc.data[g][1] + "," + qc.data[g][2] + ")\n";
			}
			else if (qc.data[g][0] == "h") {
				qiskit_py += "qc.h(" + qc.data[g][1] + ")\n";
			}
			else if (qc.data[g][0] == "cx") {
				qiskit_py += "qc.cx(" + qc.data[g][1] + "," + qc.data[g][2] + ")\n";
			}
			else if (qc.data[g][0] == "m") {
				qiskit_py += "qc.measure(" + qc.data[g][1] + "," + qc.data[g][2] + ")\n";
			}
		}

		return qiskit_py;

	}

};


/*
Begin added functionality specific to EntangledStates game
Written by Nick Singstock and Elies Gil Fuster
*/

vector<vector<double>> polar_statevector(vector<std::complex<double>> cartesian) {
	// polar formula for complex numbers from cartesian
	// it translates an entire statevector, so 2^N entries
	// each of them will be a pair (magnitude, argument)
	// declare variables
	float tol_gs = .01;
	int n_gs = cartesian.size();
	double mag_gs, arg_gs, re_gs, im_gs; // auxiliary for each entry
	vector<vector<double>> polar_statevector;

	// loop through state vector
	for (int i = 0; i < n_gs; i++) {
		// get faulty results
		re_gs = real(cartesian[i]);
		im_gs = imag(cartesian[i]);
		mag = sqrt(re_gs * re_gs + im_gs * im_gs);
		if (abs(re_gs) < tol_gs) {
			if (abs(im_gs) < tol_gs) {
				arg_gs = 0;
			}
			else if (im_gs > tol_gs) {
				arg_gs = M_PI / 2;
			}
			else {
				arg_gs = -M_PI / 2;
			}
		}
		else if (re_gs > tol_gs) {
			arg_gs = atan(im_gs / re_gs);
		}
		else {
			arg = M_PI - atan(im_gs / re_gs);
		}

		// recover the ideal solutions
		if (abs(mag_gs - 1) < tol_gs) {
			mag_gs = 1;
		}
		if (abs(sqrt(2) * mag_gs - 1) < tol_gs) {
			mag_gs = 1 / sqrt(2);
		}
		if (abs(2 * mag_gs - 1) < tol_gs) {
			mag_gs = 0.5;
		}
		if (abs(2 * sqrt(2) * mag_gs - 1) < tol_gs) {
			mag_gs = 0.5 / sqrt(2);
		}

		if (abs(arg_gs) < tol_gs) {
			arg_gs = 0;
		}
		if (abs(arg_gs - M_PI) < tol_gs) {
			arg_gs = M_PI;
		}
		vector<double> aux_gs;
		aux_gs.push_back(mag_gs);
		aux_gs.push_back(arg_gs);

		// add entry to polar state vector
		polar_statevector.push_back(aux_gs);
	}

	return polar_statevector;
}

vector<vector<vector<double>>> start_level(int level_gs) {
	// gives initial state and goal state for each level
	// returns states, with states[0] initial state and state[1] goal state

	// initialize variables
	vector<vector<vector<double>>> states_gs;
	vector<vector<double>> initial_gs, goal_gs;
	double aux_gs = 1 / sqrt(2);

	switch (level_gs)	// Initial and goal states are hard coded for each level
	{					// This ain't pretty
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
			initial_gs.push_back(polar1);
			goal_gs.push_back(polar2);
		}
		break;

	case 1:
		for (int i = 0; i < 8; i++) {
			vector<double> polar1, polar2;
			polar1.push_back(0.5 * aux_gs);
			polar1.push_back(0.0);

			if (i == 6) {
				polar2.push_back(aux_gs);
				polar2.push_back(0.0);
			}
			else if (i == 7) {
				polar2.push_back(aux_gs);
				polar2.push_back(M_PI);
			}
			else {
				polar2.push_back(0.0);
				polar2.push_back(0.0);
			}
			initial_gs.push_back(polar1);
			goal_gs.push_back(polar2);
		}
		break;

	case 2:
		for (int i = 0; i < 8; i++) {
			vector<double> polar1, polar2;
			if (i == 0 || i == 3) {
				polar1.push_back(aux_gs);
			}
			else {
				polar1.push_back(0.0);
			}
			polar1.push_back(0.0);

			if (i == 3) {
				polar2.push_back(aux_gs);
				polar2.push_back(-0.5 * M_PI);
			}
			else if (i == 5) {
				polar2.push_back(aux_gs);
				polar2.push_back(0.5 * M_PI);
			}
			else {
				polar2.push_back(0.0);
				polar2.push_back(0.0);
			}
			initial_gs.push_back(polar1);
			goal_gs.push_back(polar2);
		}
		break;
	}

	states_gs.push_back(initial_gs);
	states_gs.push_back(goal_gs);

	return states_gs;
}

vector<string> get_secret_gates(int level_gs) {
	// let the game manager know how to generate the initial state
	// again, this is hard coded so not very beautiful

	vector<string> gates_gs;

	switch (level_gs) {
	case 1:
		gates_gs.push_back("h0");
		gates_gs.push_back("h1");
		gates_gs.push_back("h2");
		break;

	case 2:
		gates_gs.push_back("h1");
		gates_gs.push_back("c10");
		break;

	}

	return gates_gs;
}

vector<std::complex<double>> get_end_state(int level_gs) {
	// let the game manager know what the goal state is
	// the end state for each level is hard coded
	vector<std::complex<double>> end_state_gs;
	double aux_gs = 1 / sqrt(2);

	switch (level_gs) {
	case 0:
		for (int i = 0; i < 8; i++) {
			if (i == 6) {
				std::complex<double> e(-1.0, 0.0);
				end_state_gs.push_back(e);
			}
			else {
				std::complex<double> e(0.0, 0.0);
				end_state_gs.push_back(e);
			}
		}
		break;

	case 1:
		for (int i = 0; i < 8; i++) {
			if (i == 6) {
				std::complex<double> e(aux, 0.0);
				end_state_gs.push_back(e);
			}
			else if (i == 7) {
				std::complex<double> e(-aux, 0.0);
				end_state_gs.push_back(e);
			}
			else {
				std::complex<double> e(0.0, 0.0);
				end_state_gs.push_back(e);
			}
		}
		break;

	case 2:
		for (int i = 0; i < 8; i++) {
			if (i == 3) {
				std::complex<double> e(0.0, -aux);
				end_state_gs.push_back(e);
			}
			else if (i == 5) {
				std::complex<double> e(0.0, aux);
				end_state_gs.push_back(e);
			}
			else {
				std::complex<double> e(0.0, 0.0);
				end_state_gs.push_back(e);
			}
		}
		break;
	}

	return end_state_gs;
}

//int mypow(int x, int p)
//{
//	if (p == 0) return 1;
//	if (p == 1) return x;
//
//	int tmp = mypow(x, p / 2);
//	if (p % 2 == 0) return tmp * tmp;
//	else return x * tmp * tmp;
//}
//
//vector<int> x_permut(int n, int m) {
//	// performs permutation equivalent to X gate
//	// args:
//	//	n (int): number of qubits
//	//	m (int): qubit on which gate is applied
//	// rets:
//	//	perm (vector int): permutation of 2**n elements
//
//	// initialize variables
//	vector<int> perm;
//	for (int i = 0; i < mypow(2, n); i++) {
//		perm.push_back(i);
//	}
//
//	// loop through entries
//	// avoid matrix multiplication, quite fancy footwork
//	for (int i = 0; i < mypow(2, n - 1); i++) {
//		int I = 2 * i - i % mypow(2, m);
//		int J = I + mypow(2, m);
//		int a = perm[I];
//		perm[I] = perm[J];
//		perm[J] = a;
//	}
//
//	return perm;
//}
//
//vector<int> cx_permut(int n, int c, int t) {
//	// performs permutation equivalent to CX gate
//	// args:
//	//	n (int): number of qubits
//	//	c (int): control qubit
//	//	t (int): target qubit
//	// rets:
//	//	perm (vector int): permutation of 2**n elements
//	vector<int> perm;
//	for (int i = 0; i < mypow(2, n); i++) {
//		perm.push_back(i);
//	}
//
//	// cf. x_permut
//	for (int i = 0; i < mypow(2, n - 2); i++) {
//		int I = mypow(2, c) + i % mypow(2, c) + ((i - i % mypow(2, c)) * 2) % mypow(2, t) + 2 * ((i - i % mypow(2, c)) * 2 - ((2 * (i - i % mypow(2, c)) % mypow(2, t))));
//		int J = I + mypow(2, t);
//		int a = perm[I];
//		perm[I] = perm[J];
//		perm[J] = a;
//	}
//
//	return perm;
//}


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


// declare globals for game_manager:
int level_gm;
vector<string> gate_list_gm;
vector<string> secret_gates_gm;
vector<complex<double>> end_state_gm;
vector<string> full_gate_list_gm;
vector<complex<double>> current_state_gm;

// will manage the current status of the game
class game_manager {
public:
	int n_qubits = 3; // hard coded for now but could come from level

	int level_gm = 1;
	vector<string> gate_list_gm;

	// get secret states and end state goal, as defined by the level
	vector<string> secret_gates_gm = get_secret_gates(level_gm);
	vector<std::complex<double>> end_state_gm = get_end_state(level_gm);

	// combine secret_gates and gate_list to get full_gate_list
	vector<string> full_gate_list = combine_gates(secret_gates_gm, gate_list_gm);

	// get current statevector
	vector<std::complex<double>> current_state_gm = quick_states(n_qubits, full_gate_list_gm);

	// call to determine if game is won
	bool game_won() {
		// get current statevector
		current_state_gm = quick_states(n_qubits, full_gate_list_gm);
		return check_states(current_state_gm, end_state_gm);
	}

	void print_level() {
		cout << "Current level: " << level_gm << endl;
	}

	void update_level(int new_level) {
		level_gm = new_level;
		// get secret states and end state goal, as defined by the level
		secret_gates_gm = get_secret_gates(level_gm);
		end_end_state_gmstate = get_end_state(level_gm);
		// combine secret_gates and gate_list to get full_gate_list
		full_gate_list_gm = combine_gates(secret_gates_gm, gate_list_gm);
		// get current statevector
		current_state_gm = quick_states(n_qubits, full_gate_list_gm);
	}

	void add_gate(string new_gate) {
		gate_list_gm.push_back(new_gate);
		// combine secret_gates and gate_list to get full_gate_list
		full_gate_list_gm = combine_gates(secret_gates_gm, gate_list_gm);
		// get current statevector
		current_state_gm = quick_states(n_qubits, full_gate_list_gm);
	}

	void update_gates(vector<string> new_gate_list) {
		gate_list_gm = new_gate_list;
		// combine secret_gates and gate_list to get full_gate_list
		full_gate_list_gm = combine_gates(secret_gates_gm, gate_list_gm);
		// get current statevector
		current_state_gm = quick_states(n_qubits, full_gate_list_gm);
	}

	// call to get current statevector of game
	vector<vector<double>> statevector() {
		return polar_statevector(current_state_gm);
	}

	// call to get goal state
	vector<vector<double>> goal_state() {
		return polar_statevector(end_state_gm);
	}
	// call to get complex goal state for testing purposes
	vector<std::complex<double>> goal_state_complex() {
		return end_state_gm;
	}

};


UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class MYPROJECT9_API UMicroQiskitWrapperActorComponent : public UActorComponent
{
	GENERATED_BODY()

public:	
	// Sets default values for this component's properties
	UMicroQiskitWrapperActorComponent();

protected:
	// Called when the game starts
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) override;
	//UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments")
	//	TMap<FIntPoint, FVector2D> GetStateKet(int InputOperationIndex);
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<FVector2D> GetStateKet(int InputOperationIndex);
	UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<int> x_permut(int n, int m);
    UFUNCTION(BlueprintCallable, Category = "FxnGetMultipoleMoments") TArray<int> cx_permut(int n, int c, int t);
};
