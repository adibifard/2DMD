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
		Na = 1000;
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
	double T = 300; // Kelvin
	
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

	MPI_Datatype MPItwoDVecDouble;
	MPI_Type_contiguous(2, MPI_DOUBLE, &MPItwoDVecDouble);
	MPI_Type_commit(&MPItwoDVecDouble);

	MPI_Datatype MPItwoDVecInt;
	MPI_Type_contiguous(2, MPI_INT, &MPItwoDVecInt);
	MPI_Type_commit(&MPItwoDVecInt);

	MPI_Datatype MPIVectInt;
	MPI_Type_contiguous(Na, MPI_INT, &MPIVectInt);
	MPI_Type_commit(&MPIVectInt);
	
	// MPI data type for atom-type structure
	int const count = 9; // number of members in the structure
	int blocklengths[count] = { 1, 2 ,2,2,2,2,2,1,1}; // length of each block of the structure
	MPI_Aint offsets[count] = { offsetof(atom,m),offsetof(atom,pos),offsetof(atom,v),offsetof(atom,a0),offsetof(atom,a1),
		offsetof(atom,f),offsetof(atom,binIJ),offsetof(atom,w),offsetof(atom,ind)};
	MPI_Datatype AtomType[count] = { MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE ,MPI_DOUBLE ,MPI_DOUBLE ,MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT };
	MPI_Datatype MPIAtom;
	MPI_Type_create_struct(count, blocklengths, offsets, AtomType, &MPIAtom);
	MPI_Type_commit(&MPIAtom);// commit the defined MPI data type

	
	MPI_Request req;
	MPI_Status status;
	////////////////////////////////// Spatial decomposition ////////////////////////////////////////////////////////////
	

	// Ndx*Ndy <= Nproc
	int Ndx = std::ceil((double)sqrt(Nproc));  // number of spatial decomposition in the x-direction
	int Ndy = std::floor((double)Nproc / Ndx); // number of spatial decomposition in the y-direction
	double DX = LBox / Ndx;
	double DY = LBox / Ndy;
	
	std::vector<std::vector<atom>> Partitions(Nproc);
	std::vector<std::vector<atom>> GhostPart(Nproc);
	std::vector <std::vector<int>> GlobalToLocalIndex(Nproc);
	std::vector <std::vector<int>> GlobalToLocalIndexGhost(Nproc);


	// Assign memory for local atoms
	std::vector<atom> Atoms_local;
	std::vector<atom> GhostAtoms_local;
	std::vector<int>  GlobalToLocalindex_local;
	std::vector<int>  GlobalToLocalIndexGhost_local;
	std::vector<std::vector<std::vector<atom>>> Bin_local(Nb); // initialize the Bin structure
	std::vector<std::vector<std::vector<int>>> BinGhost_local(Nb); // initialize the Bin structure for the ghost atoms
	std::vector<std::vector<atom>> NeighborList_local;     // Neighboring list that comprises local atoms
	std::vector<std::vector<int>> NeighborListGhost_local; // Neighboring list that only comprises ghost atoms

	for (int i = 0; i < Nb; i++)
	{
		Bin_local[i].resize(Nb);
		BinGhost_local[i].resize(Nb);
	}
	
	if (rank == Proot)
	{
		// Initialize Atoms' velocities
		InitAtomsPos(Atoms, LBox, Na);
		InitAtomsVel(Atoms, T, Na);
		double Mom0 = SumMomentum(Atoms) * Na / Navg;

		WriteToExcel("Vx_t0.csv", Atoms, "Vx(A/fs)", 0);
		WriteToExcel("Vy_t0.csv", Atoms, "Vy(A/fs)", -1);
		WriteToExcel("x0.csv", Atoms, "x(A)", 1);
		WriteToExcel("y0.csv", Atoms, "y(A)", 2);

		// Distribute the atoms between the MPI processes
		SpatialDecomp(Atoms, Nproc, DX, DY, Ndx, Ndy, GlobalToLocalIndex,
			Partitions, GhostPart, GlobalToLocalIndexGhost);
		

	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank==Proot)
	{
		int DataSize;
		Atoms_local.insert(Atoms_local.end(), Partitions[0].begin(), Partitions[0].end()); // Add owned atoms
		Atoms_local.insert(Atoms_local.end(), GhostPart[0].begin(), GhostPart[0].end());   // Add the ghost atoms
		NeighborList_local.resize(Partitions[0].size()); // the Neighboring list should only contain the owned atoms

		GlobalToLocalindex_local = GlobalToLocalIndex[0];
		GhostAtoms_local = GhostPart[0];
		GlobalToLocalIndexGhost_local = GlobalToLocalIndexGhost[0];

		for (int r=1;r<Nproc;r++)
		{
			// Send the partitioned atomic data
			DataSize = Partitions[r].size();
			MPI_Send(&DataSize, 1, MPI_INT, r, r+1, MPI_COMM_WORLD);// send the size of the data
			MPI_Send(Partitions[r].data(), DataSize, MPIAtom, r, r+2, MPI_COMM_WORLD); // send the data itself

			// Send the global to local index information (main atoms)
			DataSize = GlobalToLocalIndex[r].size();
			MPI_Send(&DataSize, 1, MPI_INT, r, r+3, MPI_COMM_WORLD);// send the size of the data
			MPI_Send(GlobalToLocalIndex[r].data(), DataSize, MPI_INT, r, r+4, MPI_COMM_WORLD);

			// Send the ghost atoms
			DataSize = GhostPart[r].size();
			MPI_Send(&DataSize, 1, MPI_INT, r, r + 5, MPI_COMM_WORLD);// send the size of the data
			MPI_Send(GhostPart[r].data(), DataSize, MPIAtom, r, r + 6, MPI_COMM_WORLD);

			// Send the global to local index information (ghost atoms)
			DataSize = GlobalToLocalIndexGhost[r].size();
			MPI_Send(&DataSize, 1, MPI_INT, r, r + 7, MPI_COMM_WORLD);// send the size of the data
			MPI_Send(GlobalToLocalIndexGhost[r].data(), DataSize, MPI_INT, r, r + 8, MPI_COMM_WORLD);

			//std::cout << "sizeGlobal: " << Partitions[r].size() << "\n";
		}
	
	}
	else if (rank!=Proot)
	{
		// Receive the partitioned atomic data
		int DataSize_local;
		MPI_Recv(&DataSize_local, 1, MPI_INT, Proot, rank+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the size of the data
		Atoms_local.resize(DataSize_local);
		MPI_Recv(Atoms_local.data(), DataSize_local, MPIAtom, Proot, rank+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the data itself
		NeighborList_local.resize(Atoms_local.size()); // Resize the neighbronig list (only based on the owned atoms per process)

		// Receive the global to local index information (main atoms)
		MPI_Recv(&DataSize_local, 1, MPI_INT, Proot, rank+3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the size of the data
		GlobalToLocalindex_local.resize(DataSize_local);
		MPI_Recv(GlobalToLocalindex_local.data(), DataSize_local, MPI_INT, Proot, rank + 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the data itself

		// Receive the ghost atoms
		MPI_Recv(&DataSize_local, 1, MPI_INT, Proot, rank + 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the size of the data
		GhostAtoms_local.resize(DataSize_local);
		MPI_Recv(GhostAtoms_local.data(), DataSize_local, MPIAtom, Proot, rank + 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the data itself
		Atoms_local.insert(Atoms_local.end(), GhostAtoms_local.begin(), GhostAtoms_local.end()); // All atoms in process p=owned atoms+ghost atoms

		// Receive the global to local index information (ghost atoms)
		MPI_Recv(&DataSize_local, 1, MPI_INT, Proot, rank + 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the size of the data
		GlobalToLocalIndexGhost_local.resize(DataSize_local);
		MPI_Recv(GlobalToLocalIndexGhost_local.data(), DataSize_local, MPI_INT, Proot, rank + 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the data itself
		
		//std::cout << "sizeLocal: " << Atoms_local.size() << "\n";

		//std::cout << "localatom(end)x: " << Atoms_local[Atoms_local.size()-1].pos[1] << "\n";
	}

	MPI_Barrier(MPI_COMM_WORLD);

	std::cout << "Rank: " << rank << "\n";

	// Unit conversions
	double KcaltoJoule = 4186.8;
	double Jtorealunits = 1e-7;
	double unitConv = KcaltoJoule * Jtorealunits;
	double rKineticEnergy;
	double rKE_KcalPerMole;
	
	// loop over the number of time steps
	for (size_t i = 0; i <= NTsteps; i++)
	{
		
		if (i==0)
		{
			// Divide particles into bins (By Meisam)
			BinParticles(Atoms_local, binSize, Bin_local); // for owned atoms
			// Find the neighboring list
			Neighboring(Atoms_local, Bin_local, GlobalToLocalindex_local, NeighborList_local); // for owned atoms (the neighboring can contain ghost atoms as well)
			// Apply pair-wise forces (By Meisam)
			ApplyForce(Atoms_local, GlobalToLocalindex_local, NeighborList_local, eps, sigma);
		}

		double KineticEnergy = 0;

		//  Compute accelerations from forces at current position (do calculations locally)
		for (atom& atom1 : Atoms_local)
		{
			// compute a(t)
			atom1.a0 = (atom1.f / atom1.m) * unitConv;
			// compute x(t+dt)
			atom1.pos = atom1.pos + atom1.v * dt + atom1.a0 * pow(dt, 2) * 0.5;
			// Apply Periodic Boundary Condition over the position of the particles leaving the box boundaries
			PBC(atom1.pos[0]);
			PBC(atom1.pos[1]);
		}

		// Remove the current bins
		for (int i = 0; i < Bin_local.size(); i++)
		{
			for (int j = 0; j < Bin_local[0].size(); j++)
			{
				Bin_local[i][j].clear();
			}
		}

		// Update the bins
		BinParticles(Atoms_local, binSize, Bin_local);

		// Find the neighboring atoms for the atoms owned by the process of rank=rank
		Neighboring(Atoms_local, Bin_local, GlobalToLocalindex_local, NeighborList_local);

		// Update force at t+dt: f(t+dt)
		ApplyForce(Atoms_local, GlobalToLocalindex_local, NeighborList_local, eps, sigma);
		
		// Determine a(t+dt) and v(t+dt)
		for (atom& atom1 : Atoms_local)
		{
			// compute a(t+dt)
			atom1.a1 = (atom1.f / atom1.m) * unitConv;
			// compute v(t+dt)
			atom1.v = atom1.v + (atom1.a0 + atom1.a1) * dt * 0.5;
		}


		CalcInstanKE(Atoms_local, KineticEnergy);
		double KE_KcalPerMole = KineticEnergy * 1e7 / KcaltoJoule;
		//std::cout << "i: " << i << "; KE_local (Kcal/mole): " << KE_KcalPerMole << "\n";
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Reduce(&KE_KcalPerMole, &rKE_KcalPerMole, 1, MPI_DOUBLE, MPI_SUM, Proot, MPI_COMM_WORLD);

		//VelVerlet(Atoms, dt, eps, sigma,Bin, binSize);
		//std::cout << "NAtoms in rank " <<rank<<"= " << Atoms_local.size() << "\n";
		

		////////////////////////// Re-partition the atoms //////////////////////////
		if (rank==Proot)
		{
			
			std::cout << "i: " << i << "; KE (Kcal/mole): " << rKE_KcalPerMole << "\n";
			// Gather all atomic data from processes

			for (int r=0;r<Nproc;r++)
			{
				Partitions[r].clear();
				GlobalToLocalIndex[r].clear();
				GhostPart[r].clear();
				GlobalToLocalIndexGhost[r].clear();
			}

			// Distribute the atoms between the MPI processes
			SpatialDecomp(Atoms, Nproc, DX, DY, Ndx, Ndy, GlobalToLocalIndex,
				Partitions, GhostPart, GlobalToLocalIndexGhost);
		}

		// Clear all the local data
		Atoms_local.clear();
		GlobalToLocalindex_local.clear();
		GhostAtoms_local.clear();
		GlobalToLocalIndexGhost_local.clear();
		NeighborList_local.clear();

		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == Proot)
		{
			int DataSize;
			Atoms_local.insert(Atoms_local.end(), Partitions[0].begin(), Partitions[0].end()); // Add owned atoms
			Atoms_local.insert(Atoms_local.end(), GhostPart[0].begin(), GhostPart[0].end());   // Add the ghost atoms
			NeighborList_local.resize(Partitions[0].size()); // the Neighboring list should only contain the owned atoms

			GlobalToLocalindex_local = GlobalToLocalIndex[0];
			GhostAtoms_local = GhostPart[0];
			GlobalToLocalIndexGhost_local = GlobalToLocalIndexGhost[0];

			for (int r = 1; r < Nproc; r++)
			{
				// Send the partitioned atomic data
				DataSize = Partitions[r].size();
				MPI_Send(&DataSize, 1, MPI_INT, r, r + 1, MPI_COMM_WORLD);// send the size of the data
				MPI_Send(Partitions[r].data(), DataSize, MPIAtom, r, r + 2, MPI_COMM_WORLD); // send the data itself

				// Send the global to local index information (main atoms)
				DataSize = GlobalToLocalIndex[r].size();
				MPI_Send(&DataSize, 1, MPI_INT, r, r + 3, MPI_COMM_WORLD);// send the size of the data
				MPI_Send(GlobalToLocalIndex[r].data(), DataSize, MPI_INT, r, r + 4, MPI_COMM_WORLD);

				// Send the ghost atoms
				DataSize = GhostPart[r].size();
				MPI_Send(&DataSize, 1, MPI_INT, r, r + 5, MPI_COMM_WORLD);// send the size of the data
				MPI_Send(GhostPart[r].data(), DataSize, MPIAtom, r, r + 6, MPI_COMM_WORLD);

				// Send the global to local index information (ghost atoms)
				DataSize = GlobalToLocalIndexGhost[r].size();
				MPI_Send(&DataSize, 1, MPI_INT, r, r + 7, MPI_COMM_WORLD);// send the size of the data
				MPI_Send(GlobalToLocalIndexGhost[r].data(), DataSize, MPI_INT, r, r + 8, MPI_COMM_WORLD);

				//std::cout << "sizeGlobal: " << Partitions[r].size() << "\n";
			}
		}
		else if (rank != Proot)
		{
			// Receive the partitioned atomic data
			int DataSize_local;
			MPI_Recv(&DataSize_local, 1, MPI_INT, Proot, rank + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the size of the data
			Atoms_local.resize(DataSize_local);
			MPI_Recv(Atoms_local.data(), DataSize_local, MPIAtom, Proot, rank + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the data itself
			NeighborList_local.resize(Atoms_local.size());

			// Receive the global to local index information (main atoms)
			MPI_Recv(&DataSize_local, 1, MPI_INT, Proot, rank + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the size of the data
			GlobalToLocalindex_local.resize(DataSize_local);
			MPI_Recv(GlobalToLocalindex_local.data(), DataSize_local, MPI_INT, Proot, rank + 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the data itself

			// Receive the ghost atoms
			MPI_Recv(&DataSize_local, 1, MPI_INT, Proot, rank + 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the size of the data
			GhostAtoms_local.resize(DataSize_local);
			MPI_Recv(GhostAtoms_local.data(), DataSize_local, MPIAtom, Proot, rank + 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the data itself
			Atoms_local.insert(Atoms_local.end(), GhostAtoms_local.begin(), GhostAtoms_local.end()); // All atoms in process p=owned atoms+ghost atoms

			// Receive the global to local index information (ghost atoms)
			MPI_Recv(&DataSize_local, 1, MPI_INT, Proot, rank + 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the size of the data
			GlobalToLocalIndexGhost_local.resize(DataSize_local);
			MPI_Recv(GlobalToLocalIndexGhost_local.data(), DataSize_local, MPI_INT, Proot, rank + 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive the data itself

			//std::cout << "sizeLocal: " << Atoms_local.size() << "\n";

			//std::cout << "localatom(end)x: " << Atoms_local[Atoms_local.size() - 1].pos[1] << "\n";
		}

		MPI_Barrier(MPI_COMM_WORLD);

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


		//MPI_Barrier(MPI_COMM_WORLD);

		

		
	}

	MPI_Finalize();
	return 0;
}


