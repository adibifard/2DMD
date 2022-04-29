#include "MDfunctions.h"
#include <random>
#include <math.h>
#include <stdlib.h>

using std::pow;
// This is the cpp file used to code the essential functions needed for our 2DMD simulator



// A function to initialize the atoms' positions
void InitAtomsPos(std::vector<atom>& Atoms, double LBox, int Na)
{
// Uniform distribution of the atoms
//	for (size_t i=0;i<Atoms.size();i++)
//	{
//

	int sx = (int)ceil(sqrt((double)Na));
	int sy = (Na + sx - 1) / sx;

	int* shuffle = (int*)malloc(Na * sizeof(int));
	for (int i = 0; i < Na; i++)
		shuffle[i] = i;

	for (int i = 0; i < Na; i++)
	{
		//
		//  particles are not spatially sorted
		//
		long int lrand48(void);
		int j = lrand48()%(Na - i);
		int k = shuffle[j];
		shuffle[j] = shuffle[Na - i - 1];

		//
		//  particles distributed uniformly
		//
		Atoms[i].pos = {LBox * (1. + (k % sx)) / (1 + sx), LBox * (1. + (k / sx)) / (1 + sy)};

	}
	free(shuffle);
}

// A function to initialize the atoms' velocities using Maxwell-Boltzmann distribution at the set-point temperature
void InitAtomsVel(std::vector<atom> &Atoms, double T, int Na)
{
	std::random_device myseed;
	std::mt19937 Myengine(myseed());
	std::vector<double> velMag;
	for (size_t i = 0; i < Na; i++)
	{
		double SIG = std::sqrt(kB * T / Atoms[i].m);
		double Mean = 0;
		std::normal_distribution<double> NormalDist(Mean, SIG);

		Atoms[i].v = {NormalDist(Myengine),NormalDist(Myengine)};

		//velMag.push_back( std::sqrt(pow(Atoms[i].v[0],2)+ pow(Atoms[i].v[1],2)));
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
void VelVerlt(atom &Atom, double dt) 
{
	//  Compute accelerations from forces at current position

	// Euler's algorithm
	Atom.pos = Atom.pos + Atom.v * dt + Atom.f * (pow(dt, 2) / (2 * Atom.m));

	Atom.v = Atom.v + Atom.f * (dt / Atom.m);
	//periodic boundary condition
	Atom.pos += Atom.v*dt +Atom.f*0.5 * dt * dt / Atom.m;
	//if (Atom.pos < 0) {
	//		Atom.pos = Atom.pos+LBox;
         //}
	//if (Atom.pos > LBox) {
	//	Atom.pos = Atom.pos-LBox;
	//	}
  
}



// The thermostat function



// Calculate KE
double CalcInstanKE(std::vector<atom> Atoms)
{
	double KE = 0;
	for (int i = 0; i < Atoms.size(); i++) 
	{
		double V = Atoms[i].v.norm();
		KE += 0.5 * Atoms[i].m * pow(V, 2);
	}
	return KE;

}


// Calculate PE
double CalcInstanPE(std::vector<atom> Atoms)
{
	double PE = 0;
	for (int i = 0; i < Atoms.size(); i++)
	{
		PE += -Atoms[i].f.DotProd(Atoms[i].pos);
	}
	return  PE;
}



