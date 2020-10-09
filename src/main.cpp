#include "interpreter.h"
#include <vector>
#include <iostream>


int main() {


	// initialize class
	Interpreter intp;
	// set number of qubits
	intp.set_qubits(2);

	// make gates string
	vector<string> gates;
	gates.push_back("h0");
	gates.push_back("c01");


	// setup circuit
	Simulator results = intp.make_circuit(gates);

	vector<std::complex<double>> vectors = results.get_statevector();


	cout << "Statevectors from interpreter with [h0, c01]:" << endl;
	for (int j = 0; j < vectors.size(); j++) {
		cout << vectors[j] << endl;
	}

	// test quick_states
	vector<std::complex<double>> quick_vec = quick_states(2, gates);
	cout << "Statevectors from quick_states with [h0, c01]:" << endl;
	for (int j = 0; j < quick_vec.size(); j++) {
		cout << quick_vec[j] << endl;
	}

	// test check_states == true
	bool same_states = check_states(quick_vec, vectors);
	cout << "Same state (true) = " << same_states << endl;

	// test combine_gates 
	vector<string> gates2;
	gates2.push_back("h1");
	gates2.push_back("x0");
	vector<string> combined_gates = combine_gates(gates, gates2);
	cout << "combined gate order:" << endl;
	for (int j = 0; j < combined_gates.size(); j++) {
		cout << combined_gates[j] << endl;
	}

	// test 3 qubits
	vector<std::complex<double>> quick_vec3 = quick_states(3, combined_gates);
	cout << "Statevectors from quick_states for 3 qubits with [h0, c01, h1, x0]:" << endl;
	for (int j = 0; j < quick_vec3.size(); j++) {
		cout << quick_vec3[j] << endl;
	}

	// test check_states == false
	vector<std::complex<double>> quick_vec2 = quick_states(2, combined_gates);
	bool same_states2 = check_states(quick_vec2, vectors);
	cout << "Same state (false) = " << same_states2 << endl;

	// test game manager class
	game_manager gm;
	cout << "game manager loaded!" << endl;

	bool iswin = gm.game_won();
	cout << "did we win? " << iswin << endl;

	vector<vector<double>> current = gm.statevector();
	cout << "Current State:" << endl;
	for (int j = 0; j < current.size(); j++) {
		cout << "(" << current[j][0] << ", " << current[j][1] << ")" << endl;
	}

	vector<vector<double>> goal = gm.goal_state();
	cout << "Goal State:" << endl;
	for (int j = 0; j < goal.size(); j++) {
		cout << "(" << goal[j][0] << ", " << goal[j][1] << ")" << endl;
	}


	return 0;
} // clang++-7 -pthread -std=c++17 -o main main2.cpp
