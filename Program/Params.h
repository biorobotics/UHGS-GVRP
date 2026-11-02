/*  ---------------------------------------------------------------------- //
    Hybrid Genetic Search for Arc Routing Problems -- HGS-CARP
    Copyright (C) 2016 Thibaut VIDAL

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//  ---------------------------------------------------------------------- */

#ifndef PARAMS_H
#define PARAMS_H

#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <iostream>
#include "math.h"
#include <time.h>
#include <algorithm>
#include "Client.h"
#include "Vehicle.h"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <Eigen/Dense>

using namespace Eigen;
namespace py = pybind11;

typedef Matrix<long, Dynamic, 1> VectorXl;

using namespace std ;

// little function used to clear some arrays
template <class C> void FreeClear( C & cntr ) {
    for ( typename C::iterator it = cntr.begin(); 
              it != cntr.end(); ++it ) {
        delete * it;
    }
    cntr.clear();
}

// Pre-definition, to allow compilation with self-references
class Params ;
class Vehicle ;

// A little auxiliary structure
struct pairB 
{
	int myInt ;
	int iCour ; 
	Params * myparams ;
};

class Params
{
public:

	/* ------------------------- PROBLEM DATA -------------------------- */

	// random seed
	int seed ;

	int num_nodes ;
	int ar_maxNbNodesPerClient ;
	vector < vector < double > > cost_mat ;
	void ar_computeClientToClientDistances();

	// number of customers/services considered in the vehicle routing problem
	int nbClients ;

	// number of customers
	int nbTotalServices ;

	// number of vehicles per depot
	int nbVehiculesPerDep ;

	// number of depots
	const int nbDepots = 1;

	// list of vehicles available
	vector < Vehicle > ordreVehicules ;

	// number of vehicles available
	int nombreVehicules ;

	// array containing the information of each separate client/service
	Client * cli ;

	// travel time (was used for the CVRP) now its mainly used as an intermediate structure to compute the granular search proximity
	double ** timeCost ;

	// isCorrelated[i][j] returns true if and only if j is considered to be among the closest customers to i (granular search parameter)
	vector < vector <bool> > isCorrelated ;

	/* ----------------- METHOD PARAMETERS & OTHER DATA ---------------- */

	// Are we running a feasibility problem ?
	// In this case we would stop the search as soon as a feasible solution is found
	bool isSearchingFeasible ;

	// population size parameters
	int mu ; // Default 25
	int lambda ; // Default 40
	int el ; // Default 8

	// number of close individuals, taken into account in the distance measure (diversity management)
	int nbCountDistMeasure ; // Default 3

	// penalty coefficients (are adapted during the search)
	double penalityCapa ;
	double penalityLength ;

	// target of feasible individuals following LS
	double minValides ; // Default 0.25
	double maxValides ; // Default 0.35

	// number of close customers considered in RI (granular search)
	int granularity ; // Default 40

	// how much additional capacity consumption (multiplicator) allowed in Split
	double borne ; // Default 2

	// max size of a SeqData
	// The preprocessing effort could also be limited to O(n^4/3) using hierarchies as in Irnich 2008 (JOC)
	int sizeSD ; // Default 10

	/* ------------------------  PARSING ROUTINES  -------------------- */

	// setting the parameters of the method
	void setMethodParams () ;

	// reading a customer from the stream
	void getClient (int i, Client * myCli);

	// builds the other data structures (granular search etc...)
	void calculeStructures () ;
	
	// shuffle the lists of closest customers
	void shuffleProches () ;

  double cost_thresh;

	// constructor
	Params(unsigned int seedRNG, int nbVeh, double vc, const Ref<const VectorXd> &demand_per_cluster, bool take_first_soln, const std::vector<py::array_t<long>> &cluster_to_node_ptr, const Ref<const MatrixXd> &cost_mat, double cost_thresh);

	// destructor
	~Params(void);
};
#endif

