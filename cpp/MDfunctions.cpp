#include "MDfunctions.h"
#include <random>
#include <span>
#include <fstream>


using std::pow;
double LBox;
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
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
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
	double kBmod = (kB * kg2gr * pow(m2A, 2)) / pow(s2fs, 2);
	
	double vsum_x=0.0;
	double vsum_y=0.0;
	for (size_t i = 0; i < Na; i++)
	{
		double SIG = std::sqrt(kB * T / Atoms[i].m);
		double Mean = 0;
		std::normal_distribution<double> NormalDist(Mean, SIG);

		Atoms[i].v = {NormalDist(Myengine),NormalDist(Myengine)};
        vsum_x = vsum_x + Atoms[i].v[0];
		vsum_y = vsum_y + Atoms[i].v[1];
		//velMag.push_back( std::sqrt(pow(Atoms[i].v[0],2)+ pow(Atoms[i].v[1],2)));
	}
	for (size_t i = 0; i < Na; i++)
	{
	    Atoms[i].v[0]= Atoms[i].v[0]-vsum_x/Na;
	    Atoms[i].v[1]= Atoms[i].v[1]-vsum_x/Na;
	}
	vsum_x=0.0;
	vsum_y=0.0;
	for (size_t i = 0; i < Na; i++)
	{
	    vsum_x = vsum_x + Atoms[i].v[0] ;
	    vsum_y = vsum_y + Atoms[i].v[1] ;
	}
	
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

				double Fij = 4 * eps * (12 * pow(sig, 12) / pow(r, 13) - 6 * pow(sig, 6) / pow(r, 7));

				// Update the force components for the i-th particle in 2D
				Atoms[i].f = Atoms[i].f - nij * Fij;
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

	// Apply Periodic Boundary Condition over the position of the particles leaving the box boundaries
	PBC(Atom.pos[0]);
	PBC(Atom.pos[1]);
  
}

// Velocity Verlet
void VelVerlet(std::vector<atom>& Atoms, double dt, double eps, double sig)
{
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
	}

	// Update force at t+dt
	ApplyForce(Atoms, eps, sig);

	for (atom &atom1 : Atoms)
	{
		// compute a(t+dt)
		atom1.a1 = (atom1.f / atom1.m) * unitConv;
		// compute v(t+dt)
		atom1.v = atom1.v + (atom1.a0 + atom1.a1) * dt * 0.5;
		// Apply Periodic Boundary Condition over the position of the particles leaving the box boundaries
		PBC(atom1.pos[0]);
		PBC(atom1.pos[1]);
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
