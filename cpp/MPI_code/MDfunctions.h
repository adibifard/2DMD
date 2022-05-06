#pragma once
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>      
#include <cstring>

// Define atomic constants
#define kB 1.380649e-23 // J/k
#define Navg 6.0221408e+23
#define m2A 1e10
#define s2fs 1e15
#define kg2gr 1e3
#define kB_ru kB*1e-7 // (gr*A^2)/(fs^2*k)

#define rcut_off 10


// Define a structure to overload the 2D vector type
template<typename T>
struct TwoDvec
{
private:
	T x, y;
	int size = 2;

public:
	
	// initializer
	TwoDvec(T x = 0, T y = 0) :x(x), y(y) {}

	// The norm method
	T norm() { return std::sqrt(x * x + y * y); };

	// Overload + operator
	TwoDvec<T> operator+(const TwoDvec<T> rhs) {
		return TwoDvec<T>(x + rhs.x, y + rhs.y);
	}

	// Overload [] operator
	T& operator[](const int index) 
	{
		if (index >= size) {
			std::cout << "Array index out of bound, exiting";
			exit(0);
		}
		if (index == 0) {
			return x;
		} else if (index == 1) {
			return y;
		}
			

	}

	// Overload += operator
	TwoDvec<T> operator+=(const TwoDvec<T> rhs) {
		return TwoDvec<T>(x += rhs.x, y += rhs.y);
	}

	// Overload - operator
	TwoDvec<T> operator-(const TwoDvec<T> rhs) {
		return TwoDvec<T>(x - rhs.x, y - rhs.y);
	}

	// Overload -= operator
	TwoDvec<T> operator-=(const TwoDvec<T> rhs) {
		return TwoDvec<T>(x -= rhs.x, y -= rhs.y);
	}


	// Overload the == operator
	bool operator==(const TwoDvec<T> rhs) const {
		return x == rhs.x
			&& (y == rhs.y);
	}


	// Overload * for multiplication of vector with a scalar
	TwoDvec<T> operator*(const T rhs)
	{
		return TwoDvec<T>(x * rhs, y * rhs);
	}

	//  Dot product method for vector-vector dot product
	T DotProd(const TwoDvec<T> vec2)
	{
		return (x * vec2.x+ y * vec2.y);
	}

	// Overload / for division of vector by a scalar
	TwoDvec<T> operator/(const T rhs)
	{
		return TwoDvec<T>(x / rhs, y / rhs);
	}

	// Overload / for division of vector by a scalar
	TwoDvec<int> ceil(const double number)
	{
		return TwoDvec<int>(std::ceil(x / number), std::ceil(y / number));
	}

};


typedef struct atom
{
	double m; // mass
	TwoDvec<double> pos; // 2D position
	TwoDvec<double> v; // 2D velocity
	TwoDvec<double> a0; // 2D accelaration at current time-step
	TwoDvec<double> a1; // 2D accelaration at t+dt
	TwoDvec<double> f; // 2D force
	TwoDvec<int> binIJ; // (i,j) pairs of the bin containing the particle
	//std::vector<int> NeighbIndex;
};

class PBC_images
{

public:
	std::vector<atom> XplusImage; std::vector<atom> XminusImage; std::vector<atom> YplusImage; std::vector<atom> YminusImage;
	std::vector<atom> XYplusImage; std::vector<atom> XYminusImage; std::vector<atom> XplusYminusImage; std::vector<atom> XminusYplusImage;

	void ReturnImageBoxes(std::vector<atom> Atoms);

};

// read the command line inputs
int  read_cmd(int argc, char** argv, char* option, int value_default);

// Initialize the position and velocity, and checksum the total momentum
void InitAtomsPos(std::vector<atom>& Atoms, double LBox, int Na);
void InitAtomsVel(std::vector<atom> &Atoms, double T,int Na);
double SumMomentum(std::vector<atom> Atoms);

double setBoxSize(const double density, const int N);
void BinParticles(std::vector<atom>& Atoms, const double BinSize, std::vector<std::vector<std::vector<int>>>& Bin);
// Integration Algorithms
void Euler(atom &particle,double dt);
void VelVerlet(std::vector<atom>& Atoms, double dt, double eps, double sig, std::vector<std::vector<std::vector<int>>>& Bin, double binSize);
void PBC(double& x);

// Force calculations
void Neighboring(std::vector<atom>& Atoms, std::vector<std::vector<std::vector<int>>> Bin, std::vector<int> GlobalToLocalIndex, std::vector<std::vector<int>>& NeighborList_local);
void ApplyForce(std::vector<atom>& Atoms, std::vector<int> GlobalToLocalIndex, std::vector<std::vector<int>> NeighborList_local, double eps, double sig);
// functions for property calculations
void CalcInstanKE(std::vector<atom> Atoms, double& KE);
double CalcInstanPE(std::vector<atom> Atoms);

// In-Out Functions
void WriteToExcel(std::string filename, std::vector<atom>  data, std::string colName,int Dtype);
void AtomsToGRO(std::string filename, std::vector<std::vector<atom>> AlltimeAtoms, std::string BoxDimAngle);

int GetNumDigits(int N);

// MPI-related functions
void SpatialDecomp(std::vector<atom> Atoms, int Nproc, double DX, double DY, int Ndx, int Ndy, std::vector <std::vector<int>>& GlobalToLocalIndex,
	std::vector<std::vector<atom>>& Partitions, std::vector<std::vector<atom>>& GhostPart, std::vector <std::vector<int>>& GlobalToLocalIndexGhost);