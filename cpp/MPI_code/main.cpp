#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>    
#include <sstream>

#include "stddef.h"
#include "mpi.h"
#include "MDfunctions.h"

#define cutoff 10
#define ns2fs 1e+6
// This is the entering point of the program
// Units: distance=nm, time=fs, Energy=Kj, 
// Simulation of an NVT ensemble with Periodic Boundary Conditions

int main(int argc, char** argv)
{
	// Gte the input from the commandline
	/*char Nopt = '-n';
	char* ptr = &Nopt;
	int n = read_cmd(argc, argv, ptr, 100);*/
	int Na;
	if (argc>=2)
	{
		Na= std::stoi(argv[2]);
	}
	else 
	{ 
		Na = 100;
	}


	// Setup simulation time settings
	double dt = 2; // in fs
	double tf =2 ;// in ns
	size_t NTsteps = std::round(tf* ns2fs /dt);
	std::cout << "NTimeSteps: " << NTsteps << "\n";

	// Set forcefield (LJ) parameters (Argon-> from ATB molecular structures)
	double sigma = 3.409840044152281813;	 // in A
	double eps = 2.381475572067369983E-01;   // Kcal/mole
	double mass = 39.948000000;				// gr/mole

	// Set the number of atoms in the box
	
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


	

	////////////////////////////////// Initialize MPI ////////////////////////////////////////////////////////////////
	int Nproc, rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &Nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int Proot = 0;
	//////////////////////////////////// Map the structured atom data type for MPI processes ////////////////////////////////
	int const count = 8; // number of members in the structure
	int blocklengths[count] = { 1, 2 ,2,2,2,2,2,Na}; // length of each block of the structure
	MPI_Aint offsets[count] = { offsetof(atom,m),offsetof(atom,pos),offsetof(atom,v),offsetof(atom,a0),offsetof(atom,a1),offsetof(atom,f),offsetof(atom,binIJ),offsetof(atom,NeighbIndex)};
	MPI_Datatype AtomType[count] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE ,MPI_DOUBLE ,MPI_DOUBLE ,MPI_DOUBLE,MPI_INT,MPI_INT };
	MPI_Datatype MPIAtom;
	MPI_Type_create_struct(count, blocklengths, offsets, AtomType, &MPIAtom);
	MPI_Type_commit(&MPIAtom);// commit the defined MPI data type

	//////////////////////////////////// Map the structured Bin data type for MPI processes ////////////////////////////////
	/*std::vector<std::vector<std::vector<int>>> Bin(Nb);*/

	int const NumElBin = 1;
	int blockBin[NumElBin] = { 1 };
	//MPI_Aint OFFsetsBin

	MPI_Status status;
	MPI_Request req;
	////////////////////////////////// Spatial decomposition ////////////////////////////////////////////////////////////
	// Assign memory for local atoms
	std::vector<atom> Atomslocal(Na);

	// Ndx*Ndy <= Nproc
	int Ndx = std::ceil(sqrt(Nproc));  // number of spatial decomposition in the x-direction
	int Ndy = std::floor(Nproc / Ndx); // number of spatial decomposition in the y-direction
	double DX = LBox / Ndx;
	double DY = LBox / Ndy;

	std::vector<std::vector<atom>> Partitions(Nproc);
	std::vector<std::vector<atom>> GhostPart(Nproc);

	std::vector<int> particlesPerProc(Nproc);
	std::vector <std::vector<int>> GlobalToLocalIndex(Nproc);

	if (rank == Proot) 
	{
		// Decompose the atoms into spatial blocks
		for (int i = 0; i < Atoms.size(); i++) 
		{
			int I = std::ceil(Atoms[i].pos[0] / DX); 
			int J = std::ceil(Atoms[i].pos[1] / DY);
			I = (I == 0) ? I += 1 : I;
			J = (J == 0) ? J += 1 : J;
			int ProcRank = (J-1) * Ndx + I-1;
			//std::cout <<i<<" " << "; Dx: " << DX <<"Dy: "<<DY<< " " << I << J << " " << "ProcRank: " << ProcRank << "\n";
			GlobalToLocalIndex[ProcRank].push_back(i);
			Partitions[ProcRank].push_back(Atoms[i]);
		}
		
		for (int r = 1; r < Nproc; r++)
		{
			int DataSize = Partitions[r].size();
			//MPI_Ibcast(&Partitions[r], DataSize, MPIAtom, 0, MPI_COMM_WORLD, &req);
			MPI_Isend(&Partitions[r], DataSize, MPIAtom, r, 0, MPI_COMM_WORLD, &req);
			MPI_Wait(&req, &status);
			MPI_Isend(&GlobalToLocalIndex[r], DataSize, MPIAtom, r, 0, MPI_COMM_WORLD, &req);
			MPI_Wait(&req, &status);
			MPI_Isend(&Bin, Nb * Nb * Na, MPI_INT, r, 0, MPI_COMM_WORLD, &req);
			MPI_Wait(&req, &status);
			/*
			MPI_Send(&Partitions[r], Na, MPIAtom, r, 0, MPI_COMM_WORLD);
			MPI_Recv(&Atomslocal, Na, MPIAtom, Proot, 0, MPI_COMM_WORLD, &status);*/
		}
		
		

		//std::cout << "atoms in the first processor\n";
		//for (atom atom : Partitions[0]) std::cout << "x: " << atom.pos[0] << " y: " << atom.pos[1] << "\n";
		
		
	}


	/*int particle_per_proc = (Na + Nproc - 1) / Nproc;
	int partition_offsets[Nproc + 1];
	for (int i = 0; i < Nproc + 1; i++) partition_offsets[i] = std::min((i * particle_per_proc), Na);
	std::vector<int> partition_sizes(Nproc );
	for (int i = 0; i < Nproc ; i++) partition_sizes[i] = partition_offsets[i+1]- partition_offsets[i];

	MPI_Scatterv(&Atoms, &partition_sizes, &partition_offsets, MPIAtom, &Atomslocal, Na, MPIAtom, 0, MPI_COMM_WORLD);*/

	
	std::cout << "Rank: " << rank << "\n";

	//MPI_Bcast(&Atoms, Na, MPIAtom, 0, MPI_COMM_WORLD);

	
	// Unit conversions
	double KcaltoJoule = 4186.8;
	double Jtorealunits = 1e-7;
	double unitConv = KcaltoJoule * Jtorealunits;
	double rKineticEnergy;

	
	// loop over the number of time steps
	for (size_t i = 0; i <= NTsteps; i++)
	{
		if (i==0)
		{
			// Divide particles into bins (By Meisam)
			BinParticles(Partitions[rank], binSize, Bin);
			// Find the neighboring list
			Neighboring(Partitions[rank], Bin, GlobalToLocalIndex[rank]);
			// Apply pair-wise forces (By Meisam)
			ApplyForce(Partitions[rank], eps, sigma);
		}

		double KineticEnergy = 0;
		// Populate the image boxes
		//IMGBox.ReturnImageBoxes(Atoms);

		// retrieve the spatial i and j of each process

		int Jp = std::ceil((rank + 1) / Ndx);
		
		int Ip = (rank + 1) - (Jp - 1) * Ndx;

		std::cout << "Rank: " << rank << "Ip: " << Ip << "; Jp: " << Jp << "\n";
		std::cout << "atoms in the "<<rank<<"  processor\n";
		for (atom atom : Partitions[rank]) std::cout << "x: " << atom.pos[0] << " y: " << atom.pos[1] << "\n";
	

		//
		// Update per atomic (r,v) data locally (for each process)
		// 

		//  Compute accelerations from forces at current position (do calculations locally)
		for (atom& atom1 : Partitions[rank])
		{
			// compute a(t)
			atom1.a0 = (atom1.f / atom1.m) * unitConv;
			// compute x(t+dt)
			atom1.pos = atom1.pos + atom1.v * dt + atom1.a0 * pow(dt, 2) * 0.5;
			// Apply Periodic Boundary Condition over the position of the particles leaving the box boundaries
			PBC(atom1.pos[0]);
			PBC(atom1.pos[1]);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		// Remove the current bins
		for (int i = 0; i < Bin.size(); i++)
		{
			for (int j = 0; j < Bin[0].size(); j++)
			{
				Bin[i][j].clear();
			}
		}


		// Update bins
		if (rank == Proot) 
		{
			BinParticles(Atoms, binSize, Bin);
			for (int r = 1; r < Nproc; r++)
			{
				MPI_Isend(&Bin, Nb * Nb * Na, MPI_INT, r, 0, MPI_COMM_WORLD, &req);
				MPI_Wait(&req, &status);
			}
		}
		

		Neighboring(Partitions[rank], Bin, GlobalToLocalIndex[rank]);

		// Update force at t+dt
		ApplyForce(Partitions[rank], eps, sigma);

		for (atom& atom1 : Partitions[rank])
		{
			// compute a(t+dt)
			atom1.a1 = (atom1.f / atom1.m) * unitConv;
			// compute v(t+dt)
			atom1.v = atom1.v + (atom1.a0 + atom1.a1) * dt * 0.5;
		}

		// Check if an atom has left the sub-domain of the current process

		//VelVerlet(Atoms, dt, eps, sigma,Bin, binSize);
		std::cout << "Nrank= " << Partitions[rank].size() << "\n";
		Partitions[rank].clear();
		GlobalToLocalIndex[rank].clear();
		
		// Re-partition the atoms
		if (rank == Proot)
		{
			// Decompose the atoms into spatial blocks
			for (int i = 0; i < Atoms.size(); i++)
			{
				int I = std::ceil(Atoms[i].pos[0] / DX);
				int J = std::ceil(Atoms[i].pos[1] / DY);
				I = (I == 0) ? I += 1 : I;
				J = (J == 0) ? J += 1 : J;
				int ProcRank = (J - 1) * Ndx + I - 1;
				Partitions[ProcRank].push_back(Atoms[i]);
			}
			for (int r = 1; r < Nproc; r++)
			{
				int DataSize = Partitions[r].size();
				MPI_Isend(&Partitions[r], DataSize, MPIAtom, r, 0, MPI_COMM_WORLD, &req);
				MPI_Wait(&req, &status);
				MPI_Isend(&GlobalToLocalIndex[r], DataSize, MPIAtom, r, 0, MPI_COMM_WORLD, &req);
				MPI_Wait(&req, &status);
			}
		}
	
		
		std::cout << "Nrank= " << Partitions[rank].size() << "\n";


		// Collect all the atomic data from processors (All-to-All)
		//double Mom = SumMomentum(Atoms)*Na/Navg;

		//if ((i + 1) % 500 == 0)
		//{
		//	// Determine the instantaneous pressure and temperature of the system (By Meisam and Anishumitra)
		//	 CalcInstanKE(Partitions[rank], KineticEnergy);
		//}

		CalcInstanKE(Partitions[rank], KineticEnergy);


		// Dump (r,v) data only using the root process
		//if (rank == 0)
		//{
		//	if ((i + 1) % 500 == 0)
		//	{
		//		// Determine the instantaneous pressure and temperature of the system (By Meisam and Anishumitra)
		//		double KineticEnergy = CalcInstanKE(Atoms);
		//		double KcaltoJoule = 4186.8;
		//		double KE_KcalPerMole = KineticEnergy * 1e7 / KcaltoJoule; // Kcal/mole
		//		double PEnergy = CalcInstanPE(Atoms);
		//		double TotalEn = KineticEnergy + PEnergy;
		//		std::cout << "Time (fs): " << (i + 1) * dt << "\n";
		//		std::cout << "KE: " << KE_KcalPerMole << "\n";

		//		//std::cout << "Total Energy (Kcal/mole): " << TotalEn << "\n";
		//		//std::cout << "Time (fs): "<<(i+1)*dt<<"-PE (Kcal/mole): " << PEnergy <<"-KE (Kcal/mole): "<< KE_KcalPerMole<< "\n";
		//		double kBmod = (kB * kg2gr * pow(m2A, 2)) / pow(s2fs, 2);
		//		double Tinst = KE_KcalPerMole / (Navg * kB / 4186.8);

		//		std::string XfileName = "x_t" + std::to_string(i + 1) + ".csv";
		//		std::string YfileName = "y_t" + std::to_string(i + 1) + ".csv";
		//		WriteToExcel(XfileName, Atoms, "x(A)", 1);
		//		WriteToExcel(YfileName, Atoms, "y(A)", 2);
		//		std::string fileName = "Vx_t" + std::to_string(i + 1) + ".csv";
		//		WriteToExcel(fileName, Atoms, "Vx(A/fs)", 0);

		//		//std::cout << "Tinst: " << Tinst << "\n";

		//		//Alltime_Atoms.push_back(Atoms);

		//		/*if (i == 10) {
		//			int p = 1;
		//			AtomsToGRO("trajectory5.gro", Alltime_Atoms, BoxDimAngle);
		//		}*/
		//	}
		//}


		

		MPI_Reduce(&KineticEnergy, &rKineticEnergy, 1, MPI_DOUBLE, MPI_SUM, Proot, MPI_COMM_WORLD);
		std::cout << "i: "<<i<<"; KE: " << rKineticEnergy<<"\n";

		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}


