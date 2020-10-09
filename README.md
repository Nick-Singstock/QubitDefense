# EntangledStates
 A puzzle game where the player builds circuits to unentangle a series of secretly entangled qubits. This directory is the backend for the game developed in Unreal Engine.
 
 
## How to use master_header.h from POV of game engine:

Start by importing the main library and calling the main class, `game_manager`, who will deal with everything
~~~~
#include master_header.h

game_manager gm;
~~~~

Now the game is initialized to level 1, right now we have three levels (0, 1, and 2), you can check the documentation to change levels.
At this point, you'll want to know what the initial quantum state for this level is, as well as the goal for the player.
You can achieve that using the methods `statevector()` and `state_goal()`, for example:
~~~~
vector<vector<double>> psi = gm.statevector(); # each amplitude is a magnitude-angle pair
vector<vector<double>> goal = gm.state_goal();  # same as state_vector
~~~~

This is as far as you need to go to set up the game.

Now, the player interacts with the environment and decides to place a gate (Paulis, Hadamard or CNOT) on a certain (pair of) qubit(s).
The front-end codifies the choice as a string, the first character is a letter signaling what gate is being implemented, and another (pair of) number character(s) informs on what qubit(s) the gate is applied onto.
Examples of these would be `"x0", "h2", "c10", "c02", "y1"`.

Once the input for next move is decided, the front-end calls the `add_gate` method. In this example, we will use `"x0"`:
~~~~
gm.add_gate("x0");
# Now we retrieve the current quantum state again
psi = gm.statevector();
bool iswin = gm.game_won(); # checks whether the current state equals the goal state
~~~~
With `game_won`, the front-end knows when the player has reached their objective.
