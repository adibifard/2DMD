#pragma once
#include <vector>
#include <cmath>
#include <iostream>

// Define atomic constants
#define kB 1.380649e-23 // J/k

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
	TwoDvec<double> f; // 2D force
	TwoDvec<int> binIJ; // (i,j) pairs of the bin containing the particle
};

void InitAtomsPos(std::vector<atom> &Atoms, double L_Box, int Na);
void InitAtomsVel(std::vector<atom> &Atoms, double T,int Na);
double setBoxSize(const double density, const int N);
void BinParticles(atom& particle, double BinSize);
void VelVerlt(atom &particle,double dt);
void ApplyForce(atom& Atoms_i, atom Atoms_j, double sig, double eps);

// functions for property calculations
double CalcInstanKE(std::vector<atom> Atoms)
