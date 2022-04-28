#include "MDfunctions.h"
#include <random>

using std::pow;
// This is the cpp file used to code the essential functions needed for our 2DMD simulator



// A function to initialize the atoms' positions
void InitAtomsPos(std::vector<atom>& Atoms, double LBox, int Na)
{
	// Uniform distribution of the atoms
//	for (size_t i=0;i<Atoms.size();i++)
//	{
//		
//	}
	int n, p, i, j;
	double a;
	//
	//   // Number of atoms in each direction
	n = int(ceil(pow(Atoms.size(), 1.0 / 3)));
	//
	//   //  spacing between atoms along a given direction
	a = LBox / n;
	//   
	//   //  index for number of particles assigned positions
	p = 0;
	//   //  initialize positions
	for (i=0; i<n; i++) {
     for (j=0; j<n; j++) {
		 if (p < Na) {
			 Atoms[p].pos = { (i)*a,(j)*a }; //atom at the origin
				 p = p + 1;
			 Atoms[p].pos = { (i)*a,(j + 0.5) * a }; //atom at one edge
				 p = p + 1;
			 Atoms[p].pos = { (i + 0.5) * a,(j)*a }; //atom at the other edge
				 p = p + 1;
			 Atoms[p].pos = { (i + 0.5) * a,(j + 0.5) * a }; //atom diagonal from the origin
				 p = p + 1;
		 }
       }
     }
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
	if (Atom.pos < 0) {
			Atom.pos = Atom.pos+LBox;
         }
	if (Atom.pos > LBox) {
		Atom.pos = Atom.pos-LBox;
		}
  
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



