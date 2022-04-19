#pragma once
#include <vector>

typedef struct atom
{
	std::vector<double> r;
	std::vector<double> v;
	std::vector<double> f;
	std::vector<int> binIJ;
};

void VelVerlt(atom &particle);
void ApplyForce(atom& Atoms_i, atom& Atoms_j);
