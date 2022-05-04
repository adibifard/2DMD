#include "MDfunctions.h"
#include <random>
#include <span>
#include <fstream>
#include <sstream>
#include <string>

#include "mpi.h"

using std::pow;
double LBox;
// This is the cpp file used to code the essential functions needed for our 2DMD simulator


// read the command line inputs
int  read_cmd(int argc, char** argv, char* option, int value_default)
{
	for (int i = 0; i < argc; i++) 
	{
		if (strcmp(argv[i],option) == 0)
		{
			return std::stoi(argv[i + 1]);
		}
	}
	
	return value_default;
}


// A function to initialize the atoms' positions
void InitAtomsPos(std::vector<atom>& Atoms, double LBox, int Na)
{
	// Uniform distribution of the atoms
//	for (size_t i=0;i<Atoms.size();i++)
//	{
//		
//	}



	/*int n, p, i, j;
	double a;*/
	//
	int Nxy = ceil(sqrt(Na));
	double dxy = LBox / Nxy;

	int p = 0;
	for (int i = 0; i < Nxy; i++) {
		for (int j = 0; j < Nxy; j++) {
			if (p < Na-1) {
				Atoms[p].pos = { i * dxy , j * dxy };
				p += 1;
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
	double kBmod = (kB * kg2gr * pow(m2A, 2)) / pow(s2fs, 2);

	double vsum_x = 0.0;
	double vsum_y = 0.0;
	for (size_t i = 0; i < Na; i++)
	{
		double SIG = std::sqrt(kB_ru * T / (Atoms[i].m/ Navg));// kb_ru=[gr.A^2/(fs^2*K)]

		double Mean = 0;
		std::normal_distribution<double> NormalDist(Mean, SIG);

		Atoms[i].v = {NormalDist(Myengine),NormalDist(Myengine)};

		//velMag.push_back( std::sqrt(pow(Atoms[i].v[0],2)+ pow(Atoms[i].v[1],2)));
	}
	for (size_t i = 0; i < Na; i++)
	{
		Atoms[i].v[0] = Atoms[i].v[0] - vsum_x / Na;
		Atoms[i].v[1] = Atoms[i].v[1] - vsum_y / Na;
	}
	vsum_x = 0.0;
	vsum_y = 0.0;
	for (size_t i = 0; i < Na; i++)
	{
		vsum_x = vsum_x + Atoms[i].v[0];
		vsum_y = vsum_y + Atoms[i].v[1];
	}

}

// A function to calculate sum of momentums
double SumMomentum(std::vector<atom> Atoms) 
{
	double Smoment = 0;
	for (atom atom : Atoms) 
	{
		double V = atom.v.norm();
		Smoment += atom.m * V;
	}

	return Smoment;
}

// This function determines the size of the box for n number of atoms
double setBoxSize(const double density, const int N)
{
	LBox = std::sqrt( N/ density);
	return LBox;
}

// This function determines the bining index of each particle 
void BinParticles(std::vector<atom>& Atoms, const double BinSize, std::vector<std::vector<std::vector<int>>>& Bin)
{
	for (int i=0;i<Atoms.size(); i++)
	{
		int I = std::floor(Atoms[i].pos[0] / BinSize);
		int J = std::floor(Atoms[i].pos[1] / BinSize);
		Atoms[i].binIJ = { I,J };
		Bin[I][J].push_back(i);
	}	
	
}

void PBC_images::ReturnImageBoxes (std::vector<atom> Atoms)
{
	int Natoms = Atoms.size();
	XplusImage.resize(Natoms); XminusImage.resize(Natoms); YplusImage.resize(Natoms); YminusImage.resize(Natoms);
	XYplusImage.resize(Natoms); XYminusImage.resize(Natoms); XplusYminusImage.resize(Natoms); XminusYplusImage.resize(Natoms);

	for (size_t i = 0; i < Natoms; i++)
	{
		XplusImage[i].pos[0] = Atoms[i].pos[0] + LBox; XplusImage[i].pos[1] = Atoms[i].pos[1];
		XminusImage[i].pos[0] = Atoms[i].pos[0] - LBox; XminusImage[i].pos[1] = Atoms[i].pos[1];
		YplusImage[i].pos[0] = Atoms[i].pos[0]; YplusImage[i].pos[1] = Atoms[i].pos[1] + LBox;
		YminusImage[i].pos[0] = Atoms[i].pos[0]; YminusImage[i].pos[1] = Atoms[i].pos[1] - LBox;
		XYplusImage[i].pos[0] = Atoms[i].pos[0] + LBox; XYplusImage[i].pos[1] = Atoms[i].pos[1] + LBox;
		XYminusImage[i].pos[0] = Atoms[i].pos[0] - LBox; XYminusImage[i].pos[1] = Atoms[i].pos[1] - LBox;
		XplusYminusImage[i].pos[0] = Atoms[i].pos[0] + LBox; XplusYminusImage[i].pos[1] = Atoms[i].pos[1] - LBox;
		XminusYplusImage[i].pos[0] = Atoms[i].pos[0] - LBox; XminusYplusImage[i].pos[1] = Atoms[i].pos[1] + LBox;
	}
}

// This function determines the neighboring list per atom
void Neighboring(std::vector<atom>& Atoms, std::vector<std::vector<std::vector<int>>> Bin)
{
	
	//for (int i=0;i<Atoms.size();i++)
	//{
	//	for (int j=0;j<Atoms.size();j++)
	//	{
	//		if (Atoms[i].binIJ)
	//	}
	//
	//}

	// iterate over the bins
	for (int i=0;i<Bin.size();i++)
	{
		for (int j=0;j<Bin[0].size();j++)
		{
			for (int k = 0; k < Bin[i][j].size(); k++) 
			{
				int atomID = Bin[i][j][k];
				Atoms[atomID].NeighbIndex.insert(Atoms[atomID].NeighbIndex.end(), Bin[i][j].begin(), Bin[i][j].end()); // (i,j)

				if (i > 0) Atoms[atomID].NeighbIndex.insert(Atoms[atomID].NeighbIndex.end(), Bin[i-1][j].begin(), Bin[i-1][j].end()); // (i-1, j)
				if (i< Bin.size()-1) Atoms[atomID].NeighbIndex.insert(Atoms[atomID].NeighbIndex.end(), Bin[i+1][j].begin(), Bin[i+1][j].end()); // (i+1, j)

				if (j>0) Atoms[atomID].NeighbIndex.insert(Atoms[atomID].NeighbIndex.end(), Bin[i][j-1].begin(), Bin[i ][j-1].end()); // (i, j-1)
				if (j< Bin.size() - 1) Atoms[atomID].NeighbIndex.insert(Atoms[atomID].NeighbIndex.end(), Bin[i ][j+1].begin(), Bin[i ][j+1].end()); // (i, j+1)

				if (i < Bin.size() - 1 && j < Bin.size() - 1 ) Atoms[atomID].NeighbIndex.insert(Atoms[atomID].NeighbIndex.end(), Bin[i+1][j + 1].begin(), Bin[i+1][j + 1].end());// (i+1, j+1)
				if (i > 0 && j>0) Atoms[atomID].NeighbIndex.insert(Atoms[atomID].NeighbIndex.end(), Bin[i - 1][j-1].begin(), Bin[i - 1][j-1].end());// (i-1, j-1)

				if (i < Bin.size() - 1 && j >0) Atoms[atomID].NeighbIndex.insert(Atoms[atomID].NeighbIndex.end(), Bin[i + 1][j - 1].begin(), Bin[i + 1][j - 1].end());// (i+1, j-1)
				if (i>0 && j < Bin.size() - 1) Atoms[atomID].NeighbIndex.insert(Atoms[atomID].NeighbIndex.end(), Bin[i - 1][j + 1].begin(), Bin[i - 1][j + 1].end());// (i-1, j+1)
			}
		}
	}
}


// A function to determine the fi-- forces acting over the i-th atom
void ApplyForce(std::vector<atom>& Atoms,double eps, double sig)
{
	for (int i=0;i<Atoms.size();i++)
	{
		Atoms[i].f = { 0,0 };
		for (int j=0;j<Atoms[i].NeighbIndex.size();j++)\
		{
			int Indj = Atoms[i].NeighbIndex[j];
			// Magnitude of the force between particles i and j 
			TwoDvec<double> rij = Atoms[i].pos - Atoms[Indj].pos;

			double r = rij.norm();
			if (r != 0 && r <= rcut_off)
			{
				// The direction of the force between the particles i and j 
				TwoDvec<double> nij = rij / r;

				double Fij = 4 * eps * ((12 * pow(sig, 12) / pow(r, 13)) - (6 * pow(sig, 6) / pow(r, 7)));

				// Update the force components for the i-th particle in 2D
				Atoms[i].f -= nij * Fij;
			}
		}
	}
}


// A function to impose the PBC
void PBC(double &x)
{

	if (x < 0)
	{
		x+= LBox;
	}
	else if (x > LBox)
	{
		x-= LBox;
	}

}

// A function to impose bouncing wall condition
void BBC(atom& atom)
{

	if (atom.pos[0]<0)
	{
		atom.pos[0] = abs(atom.pos[0]);
		atom.v[0] = -atom.v[0];
	}
	else if (atom.pos[0] > LBox) {
		double dist = atom.pos[0] - LBox;
		atom.pos[0] = LBox-dist;
		atom.v[0] = -atom.v[0];
	}

	if (atom.pos[1] < 0)
	{
		atom.pos[1] = abs(atom.pos[1]);
		atom.v[1] = -atom.v[1];
	}
	else if (atom.pos[1] > LBox) {
		double dist = atom.pos[1] - LBox;
		atom.pos[1] = LBox - dist;
		atom.v[1] = -atom.v[1];
	}
}

// This function is applied over each particle separately, after force calculations.
// Basically, this function updates the current position and velocity of the particle using the calculated force per particle. 
void Euler(atom &Atom, double dt)
{
	//  Compute accelerations from forces at current position
	double KcaltoJoule = 4186.8;
	double Jtorealunits = 1e-7;
	double unitConv = KcaltoJoule * Jtorealunits;
	// Semi-implicit Euler integration
	Atom.v = Atom.v + (Atom.f / Atom.m) * dt* unitConv;
	Atom.pos = Atom.pos + Atom.v * dt;

	// Apply Periodic Boundary Condition over the position of the particles leaving the box 
	/*PBC(Atom.pos[0]);
	PBC(Atom.pos[1]);*/
	BBC(Atom);
  
}

// Velocity Verlet
void VelVerlet(std::vector<atom>& Atoms, double dt, double eps, double sig, std::vector<std::vector<std::vector<int>>>& Bin, double binSize)
{
	MPI_Init(NULL, NULL);
	int Nproc, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &Nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//  Compute accelerations from forces at current position
	double KcaltoJoule = 4186.8;
	double Jtorealunits = 1e-7;
	double unitConv = KcaltoJoule * Jtorealunits;

	for (atom &atom1 : Atoms) 
	{
		// compute a(t)
		atom1.a0 = (atom1.f / atom1.m) * unitConv;
		// compute x(t+dt)
		atom1.pos = atom1.pos + atom1.v * dt + atom1.a0 * pow(dt, 2)*0.5;
		// Apply Periodic Boundary Condition over the position of the particles leaving the box boundaries
		//BBC(atom1);
		PBC(atom1.pos[0]);
		PBC(atom1.pos[1]);
	}

	// Remove the current bins
	for (int i = 0; i < Bin.size(); i++)
	{
		for (int j = 0; j < Bin[0].size(); j++)
		{
			Bin[i][j].clear();
		}
	}

	// Update bins
	BinParticles(Atoms, binSize, Bin);
	Neighboring(Atoms, Bin);

	// Update force at t+dt
	ApplyForce(Atoms, eps, sig);

	for (atom &atom1 : Atoms)
	{
		// compute a(t+dt)
		atom1.a1 = (atom1.f / atom1.m) * unitConv;
		// compute v(t+dt)
		atom1.v = atom1.v + (atom1.a0 + atom1.a1) * dt * 0.5;
		
	}

	// Remove the current bins
	for (int i = 0; i < Bin.size(); i++)
	{
		for (int j = 0; j < Bin[0].size(); j++)
		{
			Bin[i][j].clear();
		}
	}


}



// The thermostat function


// Calculate KE
void CalcInstanKE(std::vector<atom> Atoms,double & KE)
{
	for (int i = 0; i < Atoms.size(); i++) 
	{
		double V = Atoms[i].v.norm();
		KE += 0.5 * Atoms[i].m * pow(V, 2);
	}

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


// Write data into an excel file
void WriteToExcel(std::string filename, std::vector<atom>  data, std::string colName,int Dtype)
{
	std::fstream fout;


	fout.open(filename, std::ios::out | std::ios::app);

	fout << colName << "\n";

	if (Dtype==0) // Dump vx
	{
		for (auto atom : data)
		{
			fout << atom.v[0] << "\n";
		}
	}
	if (Dtype == -1) // Dump vy
	{
		for (auto atom : data)
		{
			fout << atom.v[1] << "\n";
		}
	}
	
	if (Dtype == 1) // Dump x
	{
		for (auto atom : data)
		{
			fout << atom.pos[0] << "\n";
		}
	}

	if (Dtype == 2) // Dump y
	{
		for (auto atom : data)
		{
			fout << atom.pos[1] << "\n";
		}
	}

}

int GetNumDigits(int N)
{
	int count = 0;

	int temp = N;

	while (temp != 0) {
		count++;
		// Remove last digit of 'temp'
		temp /= 10;
	}

	return count;
}

// Write atomic data to a .gro file
void AtomsToGRO(std::string filename, std::vector<std::vector<atom>> AlltimeAtoms, std::string BoxDimAngle)
{
	std::ofstream myfile(filename);
	int decimal_places = 2;
	const double multiplier = std::pow(10.0, decimal_places);

	// Get the number of atoms
	int numAtoms = AlltimeAtoms[0].size();

	if (myfile.is_open())
	{
		for (std::vector<atom> AtomsPertime : AlltimeAtoms) 
		{
			myfile << "Atoms positions, generated by the c++ code, written by Meisam Adibifard \n";
			myfile << numAtoms << "\n";
			int numMolecCurrent = 1;
			int numAtomsCurrent = 1;

			for (atom EachAtom : AtomsPertime) 
			{
				double xO = std::ceil(EachAtom.pos[0] * multiplier/10) / multiplier;
				double yO = std::ceil(EachAtom.pos[1] * multiplier/10) / multiplier;
				double zO = std::ceil(0 * multiplier/10) / multiplier;

				int xO_numdigits = GetNumDigits(xO);
				xO_numdigits = (xO_numdigits == 0) ? xO_numdigits += 1 : xO_numdigits;
				int yO_numdigits = GetNumDigits(yO);
				yO_numdigits = (yO_numdigits == 0) ? yO_numdigits += 1 : yO_numdigits;
				int zO_numdigits = GetNumDigits(zO);
				zO_numdigits = (zO_numdigits == 0) ? zO_numdigits += 1 : zO_numdigits;

				myfile << std::fixed;
				myfile << std::setprecision(2);
				int numMolecCurrent_digits = GetNumDigits(numMolecCurrent);
				int numAtomsCurrent_digits = GetNumDigits(numAtomsCurrent);

				//// write-out atom coordinates
				for (size_t iWat = 0; iWat < (5 - numMolecCurrent_digits); iWat++)
				{
					myfile << " ";
				}
				myfile << numMolecCurrent << "SOL" << "    " << "O";

				for (size_t iOMol = 0; iOMol < (7 - numAtomsCurrent_digits); iOMol++)
				{
					myfile << " ";
				}

				myfile << numAtomsCurrent;

				for (size_t ixo = 0; ixo < (7 - xO_numdigits); ixo++)
				{
					myfile << " ";
				}
				myfile <<  xO;

				for (size_t iyo = 0; iyo < (7 - yO_numdigits); iyo++)
				{
					myfile << " ";
				}

				myfile << yO;

				for (size_t izo = 0; izo < (7 - zO_numdigits); izo++)
				{
					myfile << " ";
				}

				myfile<< zO << "\n";

				numMolecCurrent += 1;
				numAtomsCurrent += 1;
			}
			myfile << BoxDimAngle << "\n";
		}
	}
}




//////////////////////////// MPI IMPLEMENTATION FUNCTIONS /////////////////////////////////////////////
void SpatialDecomp() 
{

}