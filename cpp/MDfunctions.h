#pragma once
#include <vector>
#include <cmath>

// Define atomic constants
#define kB 1.380649e-23 // J/k

typedef struct atom
{
	double m; // mass
	std::vector<double> r; // 2D position
	std::vector<double> v; // 2D velocity
	std::vector<double> f; // 2D force
	std::vector<int> binIJ; // (i,j) pairs of the bin containing the particle
};

void InitAtomsPos(std:vector<atom>& Atoms, double L_Box);
void InitAtomsVel(std:vector<atom>& Atoms, double T);
void BinParticles(atom& particle, double BinSize);
void VelVerlt(atom &particle);
void ApplyForce(atom& Atoms_i, atom& Atoms_j);
