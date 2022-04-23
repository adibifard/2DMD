#pragma once
#include <vector>
#include <cmath>


// Define atomic constants
#define kB 1.380649e-23 // J/k

// Define a structure to overload the 2D vector type
template<typename T>
struct TwoDvec
{
	T x, y;
	T norm() { return std::sqrt(x * x + y * y); };
	// initializer
	TwoDvec(T x = 0, T y = 0) :x(x), y(y) {}

	// Overload + operator
	TwoDvec<T> operator+(const TwoDvec<T> rhs) {
		return TwoDvec<T>(x + rhs.x, y + rhs.y);
	}

	// Overload - operator
	TwoDvec<T> operator-(const TwoDvec<T> rhs) {
		return TwoDvec<T>(x - rhs.x, y - rhs.y);
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

	// Overload / for division of vector by a scalar
	TwoDvec<T> operator/(const T rhs)
	{
		return TwoDvec<T>(x / rhs, y / rhs);
	}

	// Overload / for division of vector by a scalar
	TwoDvec<int> ceil(const double number)
	{
		return TwoDvec<int>(std::ceil(x / rhs), std::ceil(y / rhs));
	}

};


typedef struct atom
{
	double m; // mass
	TwoDvec<double> pos; // 2D position
	TwoDvec<double> v; // 2D velocity
	TwoDvec<double> f; // 2D force
	TwoDvec<int> binIJ; // (i,j) pairs of the bin containing the particle
};

void InitAtomsPos(std::vector<atom> &Atoms, double L_Box);
void InitAtomsVel(std::vector<atom> &Atoms, double T);
void BinParticles(atom& particle, double BinSize);
void VelVerlt(atom &particle);
void ApplyForce(atom& Atoms_i, atom& Atoms_j);
