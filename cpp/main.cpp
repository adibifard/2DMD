#include <cmath>
#include <vector>

#include "MDfunctions.h"

#define cutoff 1

// This is the entering point of the program
// Units: distance=nm, time=fs, Energy=Kj, 
// Simulation of an NVT ensemble with Periodic Boundary Conditions
int main()
{
	// Setup simulation time settings
	double dt = 2; // in fs
	double tf =20000 ;// in fs
	size_t NTsteps = std::round(tf/dt);
	
	// Set forcefield (LJ) parameters
	double sigma =1 ; // in nm
	double eps =1 ;   // 
	double mass = 1;

	// Set the number of atoms in the box
	size_t Na = 100; // the number of atoms
	double density = 1; //  #atoms/area
	// Allocate memory for atoms
	std::vector<atom> Atoms(Na);

	for (size_t i = 0; i < Atoms.size(); i++)
		Atoms[i].m = mass;

	double LBox = setBoxSize(density, Na);

	// Set binning properties
	double alpha = 1;
	double binSize = alpha * cutoff;
	double Nb = std::ceil(LBox / binSize); // Number of bins in 1D
	int NumBins=std::pow(Nb,2); //number of total bins
	
	std::vector<std::vector<int>> Bin; // initialize the Bin structure

	// Initialize Atoms' velocities
	double T = 300; // Kelvin
	InitAtomsPos(Atoms, LBox, Na);
	InitAtomsVel(Atoms,T,Na);

	// loop over the number of time steps
	for (size_t i = 0; i <= NTsteps; i++)
	{

		// Divide particles into bins (By Meisam)
		for (size_t i = 0; i < Na; i++) 
		{
			BinParticles(Atoms[i], binSize);
		}


		// Apply pair-wise forces (By Meisam)
		for (size_t i = 0; i < Na; i++)
		{
			Atoms[i].f = { 0,0 };
			for (size_t j = 0; j < Na; j++) // this is an O(n^2) implementation that we will need to change it after implementing the Binning function into the code
			{
				ApplyForce(Atoms[i], Atoms[j],eps,sigma); //  completed
			}
		}


		// Move the particles using velocity-Verlet algorithm (By Anishumitra)
		for (size_t i = 0; i < Na; i++)
		{
			VelVerlt(Atoms[i],dt); // Please complete the vel-verlet function
		}

		// Determine the instantaneous pressure and temperature of the system (By Meisam and Anishumitra)
		double KineticEnergy = CalcInstanKE(Atoms);
		double Tinst = KineticEnergy / (Na * kB);
	}

	return 0;
}

