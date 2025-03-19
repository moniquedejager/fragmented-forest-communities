/* We simulated subcommunity dynamics in a 2-dimensional, spatial,
  semi-neutral, individual-based model. The environment consists of nX x nY
  cells, each cell having a subcommunity of n_ind individuals. Initially,
  all cells are habitable (i.e. there is no fragmentation) and, contain a
  subcommunity. Each subcommunity is populated with the highest diversity
  theoretically possible during an initialization phase and is subsequently
  allowed to stabilize, after which we run the simulations with fragmentation
  and habitat loss. For full details on the model, see De Jager et
  # al. (in prep.).

  Latest version created: 2025-03-18, by Monique de Jager */

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
void record_subcommunities();
void simulate_community_dynamics(); 
void terminateProgram();           // Stops program termination until you have pressed a key 

// global parameters:
const int sim_nr = 1;
const int nx = 45;                 // number of patches in x-direction
const int ny = 45;                 // number of patches in y-direction
const int n = 2025;                // total number of patches
const int n_ind = 1000;            // number of individuals per patch
const double mut_rate = 0.0003;    // mutation rate / probability of seedling coming from metacommunity
const bool different_dispersal = true; // whether species have different dispersal capacities, or similar ones

// parameters for habitat loss and fragmentation: 
double clustering = 2.0;           // the level of clustering of habitat patches (1 = random distribution, 5 = highly clustered)
double f_loss = 0.0;               // fraction of habitat patches lost 

// position the cells in the lattice:
vector<vector<int>> patch_type(nx, vector<int>(ny,1));            // patch type: 0 = uninhabitable matrix, 1 = subcommunity
vector<vector<vector<int>>> species(nx, vector<vector<int>>(ny, vector<int>(n_ind)));  // matrix holding the species ID-numbers per patch
vector<vector<vector<int>>> species_new(nx, vector<vector<int>>(ny, vector<int>(n_ind)));    // ID-numbers of species in the next generation

//////////////////////////////////////////////////////////////////////////////
// main program:
int main()
{
	clock_t start = clock();

	random_device rd;
	mt19937 gen(rd());
	
	// initialize patch type and species ID numbers:
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
	
	simulate_community_dynamics();

	record_subcommunities();

	clock_t end = clock();
	double elapsed = double(end - start) / CLOCKS_PER_SEC;
	cout << "Elapsed time: " << elapsed << " seconds" << endl;

	terminateProgram();
}

std::string doubleToString(double value) {
	std::ostringstream out;
	out << std::fixed << std::setprecision(2) << value; // Set fixed-point notation and precision
	return out.str();
}

void record_subcommunities()
{
	ofstream outfile;
	string filename = "./composition" + doubleToString(clustering) + "_" + doubleToString(sim_nr) + "_" + doubleToString(f_loss) + ".txt";
	outfile.open(filename, std::ios_base::app);

	for (int i = 0; i < nx; ++i) 
	{
		for (int j = 0; j < ny; ++j)
		{
			for (int k = 0; k < n_ind; ++k)
			{
				//cout << i << " ; " << j << " ; " << species[i][j][k] << endl;
				outfile << i << " ; " << j << " ; " << species[i][j][k] << endl;
			}
		}
	}
	outfile.close();
}

void simulate_community_dynamics()
{
	// simulate community dynamics: replace individuals in the subcommunities: 
	cout << "simulate community dynamics" << endl;

	std::random_device rd;
	std::mt19937 gen{ rd() };
	uniform_real_distribution<double> dist(0.0, 1.0);
	uniform_int_distribution<int> int_dist(0, n_ind - 1);

	for (int iteration = 0; iteration < 1000; ++iteration)
	{
		// every iteration, each individual (in random order) gets a chance to reproduce and disperse offspring.
		// old species ID per location will be written to species from species_new, 
		// and whether an individual has reproduced this iteration is set to false:
		vector<vector<bool>> reproduced(nx, vector<bool>(ny, false));
		for (int i = 0; i < nx; ++i)
		{
			for (int j = 0; j < ny; ++j)
			{
				for (int k = 0; k < n_ind; ++k)
				{
					species[i][j][k] = species_new[i][j][k];
				}
			}
		}

		// randomly select subcommunities and shuffle the order of individuals. Let them reproduce in that order. 
		for (int i = 0; i < n; ++i)
		{
			int ix = int_dist(gen) % nx;
			int iy = int_dist(gen) % ny;

			while (reproduced[ix][iy])
			{
				ix = (ix + 1) % nx;
				if (ix == 0) iy = (iy + 1) % ny;
			}

			// reorder individuals in subcommunity for reproduction order:
			vector<int> a = species[ix][iy];
			
			shuffle(a.begin(),a.end(), gen);

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
					species_new[new_x][new_y][r_ind] = int_dist(gen);
				}
			}
			reproduced[ix][iy] = true;
		}
		cout << "iteration number = " << iteration << endl;
	}
}

void terminateProgram()
{
	cout << "Press any character and <ENTER> to continue" << endl;
	char chAnyChar;
	cin >> chAnyChar;
	return;
}