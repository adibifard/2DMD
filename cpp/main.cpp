#include "MDfunctions.h"


// This is the entering point of the program
int main()
{

	int Na = 20; // the number of atoms

	int NTsteps = 100;
	float dt = 1; // in fs

	std::vector<atom> Atoms;

	// loop over the number of time steps
	for (size_t i = 0; i <= NTsteps; i++)
	{

		// bin the particles (By Meisam)


		// Apply pair-wise forces (By Anishumitra)
		for (size_t i = 0; i < Na; i++)
		{
			for (size_t j = 0; j < Na; j++) // this is an O(n^2) implementation that we will need to change it after implementing the Binning function into the code
			{
				ApplyForce(Atoms[i], Atoms[j]); // Please complete the vel-verlet function
			}
		}

		// Move the particles using velocity-Verlet algorithm (By Anishumitra)
		for (size_t i = 0; i < Na; i++)
		{
			VelVerlt(Atoms[i]); // Please complete the vel-verlet function
		}

		// Determine the instantaneous pressure and temperature of the system (By Meisam)


	}

	return 0;
}

