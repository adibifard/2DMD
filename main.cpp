#include <cmath>
#include <vector>
#include <string>
#include "Mdfun.h"

#define cutoff 10
#define ns2fs 1e+6
// This is the entering point of the program
// Units: distance=nm, time=fs, Energy=Kj, 
// Simulation of an NVT ensemble with Periodic Boundary Conditions
int main()
{
	// Setup simulation time settings
	double dt = 2; // in fs
	double tf = 10000 ;// in fs
	size_t NTsteps = tf;
	std::cout << "NTimeSteps: " << NTsteps << "\n";

	// Set forcefield (LJ) parameters (Argon-> from ATB molecular structures)
	double sigma = 3.409840044152281813;	 // in A
	double eps = 2.381475572067369983E-01;   // Kcal/mole
	double mass = 39.948000000;				// gr/mole

	// Set the number of atoms in the box
	size_t Na = 1000; // the number of atoms
	double density = 0.0026; //  #atoms/area ro(Ar)=1.784 gr/cm3 -> 1.784 gr/cm2
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
	for (int i = 0; i < Nb; i++)
	{
		Bin[i].resize(Nb);
	}

	// Declare a class of images
	PBC_images IMGBox;
	// Initialize Atoms' velocities
	double T = 300; // Kelvin
	InitAtomsPos(Atoms, LBox, Na);
	InitAtomsVel(Atoms,T,Na);
	double Mom0 = SumMomentum(Atoms)*Na/Navg;

	WriteToExcel("Vx_t0.csv", Atoms, "Vx(A/fs)",0);
	WriteToExcel("Vy_t0.csv", Atoms, "Vy(A/fs)", -1);
	WriteToExcel("x0.csv", Atoms, "x(A)",1);
	WriteToExcel("y0.csv", Atoms, "y(A)", 2);
	// hold all atomic data
	std::vector<std::vector<atom>> Alltime_Atoms;
	std::string BoxLenghSTR = std::to_string(LBox/10);
	std::string BoxDimAngle = BoxLenghSTR + "   " + BoxLenghSTR + "   " + BoxLenghSTR+ "   " + "0   0   0   0   0   0";
	// loop over the number of time steps
	for (size_t i = 0; i <= NTsteps; i++)
	{

		// Divide particles into bins (By Meisam)
		BinParticles(Atoms, binSize, Bin);

		// Populate the image boxes
		//IMGBox.ReturnImageBoxes(Atoms);
		

		// Apply pair-wise forces (By Meisam)

		Neighboring(Atoms, Bin);

		ApplyForce(Atoms, eps, sigma);

		VelVerlet(Atoms, dt, eps, sigma,Bin, binSize);

		double Mom = SumMomentum(Atoms)*Na/Navg;
		// Move the particles using velocity-Verlet algorithm (By Anishumitra)
		//for (size_t i = 0; i < Na; i++)
		//{
		//	Euler(Atoms[i],dt); // Please complete the vel-verlet function
		//}

		// Remove the current bins
		for (int i = 0; i < Bin.size(); i++)
		{
			for (int j=0;j<Bin[0].size();j++)
			{
				Bin[i][j].clear();
			}
		}

		/*std::string fileName="Vx_t" + std::to_string(i+1) +".csv";
		WriteToExcel(fileName, Atoms, "Vx(A/fs)",0);*/

		if ((i+1)%500==0)
		{
			std::string XfileName = "x_t" + std::to_string(i + 1) + ".csv";
			std::string YfileName = "y_t" + std::to_string(i + 1) + ".csv";
			WriteToExcel(XfileName, Atoms, "x(A)", 1);
			WriteToExcel(YfileName, Atoms, "y(A)", 2);
		
		}
		
		// Determine the instantaneous pressure and temperature of the system (By Meisam and Anishumitra)
		double KineticEnergy = CalcInstanKE(Atoms); 
		double KcaltoJoule = 4186.8;
		double KE_KcalPerMole = KineticEnergy * 1e7/KcaltoJoule; // Kcal/mole
		double PEnergy = CalcInstanPE(Atoms);
		double TotalEn = KineticEnergy + PEnergy;
		printf("Time (fs):%lf \n " ,(i + 1) * dt);
		printf("KE:%lf \n ",KE_KcalPerMole);
		
		//std::cout << "Total Energy (Kcal/mole): " << TotalEn << "\n";
		//std::cout << "Time (fs): "<<(i+1)*dt<<"-PE (Kcal/mole): " << PEnergy <<"-KE (Kcal/mole): "<< KE_KcalPerMole<< "\n";
		double kBmod = (kB * kg2gr * pow(m2A, 2)) / pow(s2fs, 2);

		double Tinst = KE_KcalPerMole / (Navg* kB/ 4186.8);

		//std::cout << "Tinst: " << Tinst << "\n";

		//Alltime_Atoms.push_back(Atoms);

		/*if (i == 10) {
			int p = 1;
			AtomsToGRO("trajectory5.gro", Alltime_Atoms, BoxDimAngle);
		}*/

	}

	return 0;
}