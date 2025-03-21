/* We simulated subcommunity dynamics in a 2-dimensional, spatial,
  semi-neutral, individual-based model. The environment consists of nX x nY
  cells, each cell having a subcommunity of n_ind individuals. Initially,
  all cells are habitable (i.e. there is no fragmentation) and, contain a
  subcommunity. Each subcommunity is populated with the highest diversity
  theoretically possible during an initialization phase and is subsequently
  allowed to stabilize, after which we run the simulations with fragmentation
  and habitat loss. For full details on the model, see De Jager et
  # al. (in prep.).

  Latest version created: 2025-03-19, by Monique de Jager */

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <windows.h>
#include <string>
#include <algorithm>
#include <vector>
#include <array>
#include <ctime>
#include <cstdlib>
#include <random>
#include <ranges>
#include <cmath> // math.h
#define PI 3.14159265;

using namespace std;

//////////////////////////////////////////////////////////////////////////////
// functions used in the model:
void initialize();
void readSimNumber();
void readParamValues();
void record_subcommunities();
void remove_habitat_patches();
void simulate_community_dynamics(); 
void terminateProgram();           // Stops program termination until you have pressed a key 

// global parameters:
const int nx = 45;                 // number of patches in x-direction
const int ny = 45;                 // number of patches in y-direction
const int n = 2025;                // total number of patches
const int n_ind = 1000;            // number of individuals per patch
const double mut_rate = 0.0003;    // mutation rate / probability of seedling coming from metacommunity
const bool different_dispersal = true; // whether species have different dispersal capacities, or similar ones
bool initialization = false; // whether this is an initialization simulation or not

// parameters for habitat loss and fragmentation: 
double clustering = 2.0;           // the level of clustering of habitat patches (1 = random distribution, 5 = highly clustered)
double f_loss = 0.5;               // fraction of habitat patches lost 
int n_patches2 = int(round(n * (1 - f_loss)));

// position the cells in the lattice:
vector<vector<int>> patch_type(nx, vector<int>(ny,1));            // patch type: 0 = uninhabitable matrix, 1 = subcommunity
vector<vector<vector<int>>> species(nx, vector<vector<int>>(ny, vector<int>(n_ind)));  // matrix holding the species ID-numbers per patch
vector<vector<vector<int>>> species_new(nx, vector<vector<int>>(ny, vector<int>(n_ind)));    // ID-numbers of species in the next generation

// simulation number:
int simNumber = 2;		// number of the row of parameter values to pick from 
int sim_nr = 1;			// number of the simulation run for a set of parameters (1 to 10)
bool simulate = true;   // whether we should continue simulations with the next set of parameters 

//////////////////////////////////////////////////////////////////////////////
// main program:
int main()
{
	while (simulate)
	{
		clock_t start = clock();

		random_device rd;
		mt19937 gen(rd());

		readSimNumber();
		readParamValues();

		if (sim_nr > 0)
		{
			initialize();
			if (initialization == FALSE) { remove_habitat_patches(); }
			simulate_community_dynamics();
			record_subcommunities();
		}
		else { simulate = false; }
		
		clock_t end = clock();
		double elapsed = double(end - start) / CLOCKS_PER_SEC;
		cout << "Elapsed time: " << elapsed << " seconds" << endl;
	}
	terminateProgram();
}

std::string doubleToString(double value) {
	std::ostringstream out;
	out << std::fixed << std::setprecision(2) << value; // Set fixed-point notation and precision
	return out.str();
}

void initialize()
{
	if (initialization)
	{
		// initialize patch type and species ID numbers of initialization simulation:
		for (int i = 0; i < nx; ++i)
		{
			for (int j = 0; j < ny; ++j)
			{
				patch_type[i][j] = 1;
				for (int k = 0; k < n_ind; ++k)
				{
					species[i][j][k] = k + 1;
					species_new[i][j][k] = k + 1;
				}
			}
		}
	}
	else {
		// initialize subcommunities using existing initialized data:

		// read in the subcommunities (x, y, species, number of individuals):
		string filename = "./community composition/composition" + doubleToString(clustering) + "_" + doubleToString(sim_nr) + "_0.00.txt";
		cout << filename << endl;
		ifstream inFile;
		inFile.open(filename, ios::in);

		if (!inFile) { // Check if the file was opened successfully
			cerr << "Error opening file." << endl;
			terminateProgram();
		}

		int x = 0;
		int y = 0;
		int spec = 0;
		int nspec = 0;
		bool reading = true;
		vector<vector<int>> id(nx, vector<int>(ny, 0));

		cout << "reading in initial subcommunities " << endl;
		while (reading) 
		{ 
			inFile >> x >> y >> spec >> nspec;
			
			//cout << x << ", " << y << ", " << spec << ", " << nspec << ", " << id[x][y] << endl;

			if (x == -1)
			{
				reading = false;
			}
			else 
			{
				for (int i = 0; i < nspec; ++i)
				{
					species[x][y][id[x][y]] = spec;
					species_new[x][y][id[x][y]] = spec;
					id[x][y] += 1;
				}
			}
		}

		inFile.close(); // Close the file

	}
	
}

void readSimNumber() 
{
	// Reads in the simulation number for validation
	ifstream inFile;
	inFile.open("inputNumber.txt", ios::in);

	inFile >> simNumber;
	inFile.close();

	simNumber = simNumber + 1;
	cout << "starting simulation number " << simNumber << endl;

	ofstream outFile;
	outFile.open("inputNumber.txt");
	outFile << simNumber;
	outFile.close();
}

void readParamValues()
{
	// Reads in the simulation number for validation
	ifstream inFile;
	inFile.open("parameter_values.txt", ios::in);

	while (simNumber > 0)
	{
		inFile >> sim_nr >> f_loss >> clustering >> initialization;
		--simNumber;
	}
	
	inFile.close();
	n_patches2 = int(round(n * (1 - f_loss)));
}

void record_subcommunities()
{
	ofstream outfile;
	string filename = "./community composition/composition" + doubleToString(clustering) + "_" + doubleToString(sim_nr) + "_" + doubleToString(f_loss) + ".txt";
	outfile.open(filename, std::ios_base::app);

	for (int i = 0; i < nx; ++i) 
	{
		for (int j = 0; j < ny; ++j)
		{
			if (patch_type[i][j] == 1)
			{
				vector<int> n_per_species(1000, 0);
				vector<int> spec = species[i][j];

				for (int k = 0; k < n_ind; ++k)
				{
					int a = spec[k] - 1;
					//cout << a << ", " << spec[k] << endl;
					n_per_species[a] += 1;
				}

				for (int k = 0; k < n_ind; ++k)
				{
					//cout << i << " ; " << j << " ; " << spec[k] << endl;
					if (n_per_species[k] > 0)
					{
						outfile << i << "	" << j << "	" << k + 1 << "	" << n_per_species[k] << endl;
					}
				}
			}
		}
	}
	// write some data to the end of the file to mark it as such:
	outfile << -1 << "	" << -1 << "	" << -1 << "	" << -1 << endl;
	outfile.close();
}

void remove_habitat_patches()
{
	// remove habitat patches, using the parameters f_loss and clustering

	// 2D Pareto distribution parameters (to calculate the weights for random sampling): 
	if (clustering == 2.0) { clustering = 2.01;  }
	double lmin = 1.00;
	double lmax = 50.0;
	double a = 0.159 * 100000.0; // equals (1/(2*pi)

	// patch types and weights:
	vector<vector<int>> h(nx, vector<int>(ny, 0));
	vector<vector<int>> w(nx, vector<int>(ny, 0));
	int n_patches = 0;
	
	std::random_device rd;
	std::mt19937 gen{ rd() };
	uniform_real_distribution<double> dist(0.0, 1.0);
	uniform_int_distribution<int> int_dist(0, nx - 1);

	// habitat is build up from one patch, until the cover level is reached:
	int ix = int_dist(gen) % nx;
	int iy = int_dist(gen) % ny;

	h[ix][iy] = 1;
	++n_patches;

	while (n_patches < n_patches2)
	{
		double sumW = 0.0;
		for (int i = 0; i < nx; ++i)
		{
			for (int j = 0; j < ny; ++j)
			{
				if ((ix == i) && (iy == j))
				{
					w[ix][iy] = 0;
				}
				else {
					double distance = sqrt(pow(ix - i, 2) + pow(iy - j, 2));
					double w_add = a * ((2 - clustering) / (pow(lmax, 2 - clustering) - pow(lmin, 2 - clustering))) * pow(distance, -1 * clustering);
					w[i][j] = w[i][j] + w_add;
				}
				sumW = sumW + w[i][j];
			}
		}
		
		// use weighted random sampling to select the next habitat patch:
		double r_value = dist(gen) * sumW;
		ix = 0;
		iy = 0;
		sumW = sumW - w[ix][iy];

		while (sumW > r_value) {
			ix += 1;
			if (ix >= nx)
			{
				ix = 0;
				iy += 1;
			}
			sumW = sumW - w[ix][iy];
		}

		h[ix][iy] = 1;
		++n_patches;

		//cout << ix << ", " << iy << ", " << n_patches << endl;
	}

	for (int i = 0; i < nx; ++i)
	{
		for (int j = 0; j < ny; ++j)
		{
			patch_type[i][j] = h[i][j];
		}
	}
}

void simulate_community_dynamics()
{
	// simulate community dynamics: replace individuals in the subcommunities: 
	cout << "simulate community dynamics" << endl;

	std::random_device rd;
	std::mt19937 gen{ rd() };
	uniform_real_distribution<double> dist(0.0, 1.0);
	uniform_int_distribution<int> int_dist(0, n_ind - 1);

	for (int iteration = 0; iteration < 10000; ++iteration)
	{
		// every iteration, each individual (in random order) gets a chance to reproduce and disperse offspring.
		// old species ID per location will be written to species from species_new, 
		// and whether an individual has reproduced this iteration is set to false:
		vector<vector<bool>> reproduced(nx, vector<bool>(ny, false));
		for (int i = 0; i < nx; ++i)
		{
			for (int j = 0; j < ny; ++j)
			{
				if (patch_type[i][j] == 1)
				{
					for (int k = 0; k < n_ind; ++k)
					{
						species[i][j][k] = species_new[i][j][k];
					}
				}
			}
		}

		// randomly select subcommunities and shuffle the order of individuals. Let them reproduce in that order.
		int n_left = n;
		while (n_left > 0)
		{
			int ix = int_dist(gen) % nx;
			int iy = int_dist(gen) % ny;

			if (patch_type[ix][iy] == 0) 
			{ 
				//cout << patch_type[ix][iy] << ", " << ix << ", " << iy << ", " << reproduced[ix][iy] << ", " << n_left << endl;
				reproduced[ix][iy] = true; 
				--n_left;
			}

			while ((reproduced[ix][iy])&&(n_left > 0))
			{
				ix = (ix + 1) % nx;
				if (ix == 0) iy = (iy + 1) % ny;

				if (patch_type[ix][iy] == 0)
				{
					//cout << patch_type[ix][iy] << ", " << ix << ", " << iy << ", " << reproduced[ix][iy] << ", " << n_left << endl;
					reproduced[ix][iy] = true;
					--n_left;
				}
			}

			if (reproduced[ix][iy] == false) {
				// reorder individuals in subcommunity for reproduction order:
				vector<int> a = species[ix][iy];

				shuffle(a.begin(), a.end(), gen);

				// reproduce and disperse: 
				for (int k = 0; k < n_ind; ++k)
				{
					double x = dist(gen);
					double lambda = static_cast<double>(a[k]) / n_ind;
					double dist_val = log10(x) / (-1 * lambda);
					double alpha = dist(gen) * 2.0 * PI;

					int new_x = (ix + static_cast<int>(round(dist_val * cos(alpha)))) % nx;
					int new_y = (iy + static_cast<int>(round(dist_val * sin(alpha)))) % ny;

					while (new_x >= nx) { new_x = nx - new_x; }
					while (new_x < 0) { new_x = nx + new_x; }
					while (new_y >= ny) { new_y = ny - new_y; }
					while (new_y < 0) { new_y = ny + new_y; }

					int r_ind = int_dist(gen);

					//cout << lambda << ", " << x << ", " << dist_val << ", " << new_x << ", " << new_y << endl;

					species_new[new_x][new_y][r_ind] = a[k];

					// there's a small chance that an individual comes in from the metacommunity:
					if (dist(gen) <= mut_rate)
					{
						species_new[new_x][new_y][r_ind] = int_dist(gen) + 1;
					}
				}
				reproduced[ix][iy] = true;
				--n_left;
			}
		}
		//cout << "iteration number = " << iteration << endl;
	}
}

void terminateProgram()
{
	cout << "Press any character and <ENTER> to continue" << endl;
	char chAnyChar;
	cin >> chAnyChar;
	return;
}