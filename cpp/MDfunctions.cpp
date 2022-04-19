#include "MDfunctions.h"


// This is the cpp file used to code the essential functions needed for our 2DMD simulator


// A function to initialize the atoms' positions
void InitAtomsPos(std:vector<atom> &Atoms, double L_Box)
{
	// Uniform distribution of the atoms

}

// A function to initialize the atoms' velocities using Maxwell-Boltzmann distribution at the specified temperature
void InitAtomsVel(std:vector<atom> &Atoms, double T)
{
	for (size_t i=0;i<Atoms.size();i++)
	{
		double alpha = std::sqrt(kB * T / Atoms[i].m);

	}
	
}

// This function determines the bining index each particle lives in
void BinParticles(atom &particle,double BinSize)
{
	particle.binIJ[0] = std::ceil(particle.r[0] / binSize);// i of the bin
	particle.binIJ[1] = std::ceil(particle.r[1] / binSize);// j of the bin
}


// This function calculates the pair-wise force between the i-th and the j-th particle
void ApplyForce(atom &Atoms_i, atom Atoms_j, double sig, double eps)
{
	// Lennard-Jone Potential
	std::vector<double> rij = Atoms_i.r - Atoms_j.r;

	Atoms_i.f += 4 * eps * (std::pow(sig/rij,12)-std::pow(sig / rij, 6));
}

// This function is applied over each particle separately, after force calculations.
// Basically, this function updates the current position and velocity of the particle using the calculated force per particle. 
void VelVerlt(atom &particle) 
{

}





// The thermostat function





