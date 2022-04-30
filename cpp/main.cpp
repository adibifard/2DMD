#include <cmath>
#include <vector>
#include <string>
#include "MDfunctions.h"

#define cutoff 10
#define ns2fs 1e+6
// This is the entering point of the program
// Units: distance=nm, time=fs, Energy=Kj, 
// Simulation of an NVT ensemble with Periodic Boundary Conditions
int main()
{
	// Setup simulation time settings
	double dt = 2; // in fs
	double tf =50 ;// in ns
	size_t NTsteps = std::round(tf* ns2fs /dt);
	
	// Set forcefield (LJ) parameters (Argon-> from ATB molecular structures)
	double sigma = 3.409840044152281813;	 // in A
	double eps = 2.381475572067369983E-01;   // Kcal/mole
	double mass = 39.948000000;				// gr/mole

	// Set the number of atoms in the box
	size_t Na = 100; // the number of atoms
	double density = 0.026; //  #atoms/area ro(Ar)=1.784 gr/cm3 -> 1.784 gr/cm2
	// Allocate memory for atoms
	std::vector<atom> Atoms(Na);
	
	
	for (size_t i = 0; i < Atoms.size(); i++)
		Atoms[i].m = mass;

	double LBox = setBoxSize(density, Na);

	// Set binning properties
	double alpha = 1;
	double binSize = alpha * cutoff;
	double Nb = std::ceil(LBox / binSize); // Number of bins in 1D
	int NumBins=std::pow(Nb,2); //number of total bins
	
	std::vector<std::vector<std::vector<int>>> Bin(Nb); // initialize the Bin structure
	for (int i = 0; i < Nb; i++) {
		Bin[i].resize(Nb);
	}
	// Declare a class of images
	PBC_images IMGBox;
	// Initialize Atoms' velocities
	double T = 350; // Kelvin
	InitAtomsPos(Atoms, LBox, Na);
	InitAtomsVel(Atoms,T,Na);
	WriteToExcel("Vx_t0.csv", Atoms, "Vx(A/fs)",0);
	//WriteToExcel("x0.csv", Atoms, "x(A)",1);
	//WriteToExcel("y0.csv", Atoms, "y(A)", 2);
	// loop over the number of time steps
	for (size_t i = 0; i <= NTsteps; i++)
	{

		// Divide particles into bins (By Meisam)
		BinParticles(Atoms, binSize, Bin);

		// Populate the image boxes
		IMGBox.ReturnImageBoxes(Atoms);
		

		// Apply pair-wise forces (By Meisam)

		Neighboring(Atoms, Bin);

		ApplyForce(Atoms, eps, sigma);

		VelVerlet(Atoms, dt, eps, sigma);

		// Move the particles using velocity-Verlet algorithm (By Anishumitra)
		//for (size_t i = 0; i < Na; i++)
		//{
		//	Euler(Atoms[i],dt); // Please complete the vel-verlet function
		//}


		// Remove current binning
		for (int i = 0; i < Bin.size(); i++)
		{
			for (int j=0;j<Bin[0].size();j++)
			{
				Bin[i][j].clear();
			}
		}

		/*std::string fileName="Vx_t" + std::to_string(i+1) +".csv";
		WriteToExcel(fileName, Atoms, "Vx(A/fs)",0);*/

		/*std::string XfileName = "x_t" + std::to_string(i + 1) + ".csv";
		std::string YfileName = "y_t" + std::to_string(i + 1) + ".csv";
		WriteToExcel(XfileName, Atoms, "x(A)", 1);
		WriteToExcel(YfileName, Atoms, "y(A)", 2);*/
		// Determine the instantaneous pressure and temperature of the system (By Meisam and Anishumitra)
		double KineticEnergy = CalcInstanKE(Atoms); 
		double KcaltoJoule = 4186.8;
		double KE_KcalPerMole = KineticEnergy * s2fs / (kg2gr * m2A* KcaltoJoule); // Kcal/mole
		double PEnergy = CalcInstanPE(Atoms);
		double TotalEn = KineticEnergy + PEnergy;
		std::cout <<i << "KE: " << KE_KcalPerMole << "\n";
		if (i == 199) {
			int p = 1;
		}
		//std::cout << "Total Energy (Kcal/mole): " << TotalEn << "\n";
		//std::cout << "Time (fs): "<<(i+1)*dt<<"-PE (Kcal/mole): " << PEnergy <<"-KE (Kcal/mole): "<< KE_KcalPerMole<< "\n";
		double kBmod = (kB * kg2gr * pow(m2A, 2)) / pow(s2fs, 2);

		double Tinst = KineticEnergy  / (Navg* kB_ru);

		std::cout << "Tinst: " << Tinst << "\n";
	}

	return 0;
}

