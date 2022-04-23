#include "MDfunctions.h"

using std::pow;
// This is the cpp file used to code the essential functions needed for our 2DMD simulator



// A function to initialize the atoms' positions
void InitAtomsPos(std::vector<atom> &Atoms, double LBox)
{
	// Uniform distribution of the atoms
	for (size_t i=0;i<Atoms.size();i++)
	{
		
	}
}

// A function to initialize the atoms' velocities using Maxwell-Boltzmann distribution at the set-point temperature
void InitAtomsVel(std::vector<atom> &Atoms, double T)
{

	for (size_t i=0;i<Atoms.size();i++)
	{
		double alpha = std::sqrt(kB * T / Atoms[i].m);

	}
	
}

// This function determines the size of the box for n number of atoms
double setBoxSize(const double density, const int N)
{
	double L = std::sqrt(density * N);
	return L;
}

// This function determines the bining index of each particle 
void BinParticles(atom &particle,const double BinSize)
{
	TwoDvec<int> tmpbinIJ = particle.pos.ceil(BinSize);
	particle.binIJ = tmpbinIJ;
	
}


// This function calculates the pair-wise force between the i-th and the j-th particle
void ApplyForce(atom& Atoms_i, atom Atoms_j, double sig, double eps)
{
	// Lennard-Jone Potential

	// Magnitude of the force between particles i and j 
	TwoDvec<double> rij = Atoms_i.pos - Atoms_j.pos;
	double r = rij.norm();
	// The direction of the force between particles i and j 
	TwoDvec<double> nij = rij / r;
	
	double Fij = 4 * eps * (12 * pow(sig, 12) / pow(r, 13) - 6 * pow(sig, 6) / pow(r, 7));

	// Update the force components for the i-th particle in 2D
	Atoms_i.f = Atoms_i.f - nij * Fij;

}

// This function is applied over each particle separately, after force calculations.
// Basically, this function updates the current position and velocity of the particle using the calculated force per particle. 
void VelVerlt(atom &Atom) 
{

}




// The thermostat function





