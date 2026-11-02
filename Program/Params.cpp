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

#include "Params.h"

void Params::setMethodParams()
{
	/* MAIN PARAMETERS OF THE METHOD */

	mu = 25 ; // Population size
	lambda = 40 ; // Number of individuals per generation
	el = 12 ; // Number of elite

	nbCountDistMeasure = 3 ; // Number of close individuals considered in the distance measure (diversity management)
	granularity = 40 ; // Restriction of the LS moves to 40 closest nodes
	minFeasibleFraction = 0.15 ; // Target proportion of feasible solution
	maxFeasibleFraction = 0.20 ; // Target proportion of feasible solution
	penalityCapa = 50 ; // Initial penalties (will evolve during the search)
	penalityLength = 50; // Initial penalties (will evolve during the search)
}

bool compPredicate(pairB i, pairB j) 
{ 
	return (i.myparams->cluster_to_cluster_cost_lb[i.myInt][i.iCour] < i.myparams->cluster_to_cluster_cost_lb[j.myInt][i.iCour]) ;
}

void Params::calculateStructures () 
{
	// Initializing other structures of the search
	vector <bool> tempB2 ;
	vector <pairB> myVector ;
	pairB myPair ;
	myPair.myparams = this ;

	// Initialization of the table of "correlation" (granular search restriction)
	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		isCorrelated.push_back(tempB2);
		for (int j=0 ; j < nbClients + nbDepots ; j++)
			isCorrelated[i].push_back(i < nbDepots || j < nbDepots);
	}

	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		cli[i].closeCustomers.clear();
		cli[i].inverseCloseCustomers.clear();
	}

	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		myVector.clear();
		for (int j=0 ; j < nbClients + nbDepots ; j++)
		{
			myPair.myInt = j ;
			myPair.iCour = i ;
			if ( i != j ) myVector.push_back(myPair);
		}

		// For each customer, sorting the list of other customers
		std::sort(myVector.begin(),myVector.end(),compPredicate);

		// And keeping only the closest ones
		for (int j=0 ; j < min(nbClients,granularity) ; j++)
		{
			cli[i].inverseCloseCustomers.push_back(myVector[j].myInt);
			cli[myVector[j].myInt].closeCustomers.push_back(i);
			isCorrelated[myVector[j].myInt][i] = true ;
		}
	}
}

void Params::getClient (int i, Client * myCli)
{
	myCli->custNum = i ;
	myCli->freq = 1 ;
	myCli->serviceDuration = 0 ;

  myCli->demand = 0 ;
  if (i == 0)
  {
    myCli->ar_nodeType = AR_DEPOT ;
    myCli->freq = 0 ;
  }
  else {
    myCli->ar_nodeType = AR_CLIENT_NODE ;
  }

	if (i >= 1) nbTotalServices += myCli->freq ;
}

// vc is vehicle capacity. Assume depot cluster is cluster 0 and depot node is node 0
Params::Params(unsigned int seedRNG, int nbVeh, double vc, const Ref<const VectorXd> &demand_per_cluster, bool take_first_soln, const std::vector<py::array_t<long>> &cluster_to_node_ptr, const Ref<const MatrixXd> &cost_mat, double cost_thresh):nbVehiculesPerDep(nbVeh), isSearchingFeasible(false), cost_thresh(cost_thresh)
{
	// Main constructor of Params
	borne = 2.0 ;
	sizeSD = 10 ;

	seed = seedRNG;
	if (seed == 0) // using the time to generate a seed when seed = 0 
		srand((unsigned int)time(NULL));
	else 
		srand(seed);
    
	// Formerly the code from preleveDonnees ();
	nbTotalServices = 0 ;

  int num_clusters = cluster_to_node_ptr.size() - 1; // Not including depot
  int num_nodes = cost_mat.rows();

	// Setting the other instance parameters
	nbClients = num_clusters;

	// Building the list of vehicles
  nombreVehicules = nbDepots*nbVehiculesPerDep;
  for (int i=0 ; i < nbDepots ; i++)
  {
    for (int j=0 ; j < nbVehiculesPerDep ; j++)
    {
      ordreVehicules.push_back(Vehicle(i,1000000,vc)); // Duration constraint set to a high value
    }
  }

	// Reading the list of customers 
	cli = new Client [nbDepots + nbClients] ;
	for (int i = 0 ; i < nbClients + nbDepots ; i++) {
		getClient (i, &cli[i]);
  }

  // Formerly the code from ar_parseOtherLinesNEARP_TP();

  num_nodes = cost_mat.rows();

	/* FIRST SPECIFY THE DEPOT CHARACTERISTICS */
	cli[0].ar_Nodes.push_back(0);
	cli[0].ar_nbNodes = 1 ;

	for (int i=1 ; i <= num_clusters; i++)
	{
    cli[i].demand = demand_per_cluster[i] ;
	}

	this->cost_mat.resize(num_nodes);
	for (int i=0 ; i < num_nodes ; i++) {
		this->cost_mat[i] = std::vector<double>(num_nodes, 1.e20);
  }

	for (int i=0 ; i < num_nodes ; i++) {
		this->cost_mat[i][i] = 0 ; 
  }

	for (int i=0 ; i < num_nodes ; i++)
	{
		for (int j=0 ; j < num_nodes ; j++)
		{
      if (i != j) {
        this->cost_mat[i][j] = cost_mat(i, j);
      }
		}
	}

	/* SPECIFYING THE NODES FOR THE CLIENT DELIVERIES */
	for (int i=1 ; i <= num_clusters ; i++)
	{
    auto ptr = cluster_to_node_ptr[i].unchecked();
		cli[i].ar_nbNodes = cluster_to_node_ptr[i].shape(0) ;
    for (int local_node_idx = 0; local_node_idx < ptr.shape(0); ++local_node_idx) {
      cli[i].ar_Nodes.push_back(ptr(local_node_idx));
    }
	}

	/* IDENTIFYING THE MAXIMUM NUMBER OF NODES FOR A DELIVERY */
	ar_maxNbNodesPerClient = 0 ;
	for (int i=0 ; i <= this->nbClients ; i++) {
		if (cli[i].ar_nbNodes > ar_maxNbNodesPerClient) {
			ar_maxNbNodesPerClient = cli[i].ar_nbNodes ;
    }
  }

  ar_computeClientToClientDistances();

	// Setting the method parameters
	setMethodParams();

	// Computing the other data structures
	calculateStructures();	
}

Params::~Params(void)
{
	for (int i=0 ; i < nbClients + nbDepots ; i++)
		delete [] cluster_to_cluster_cost_lb[i] ;
	delete [] cluster_to_cluster_cost_lb ;
	delete [] cli ;
}

void Params::shuffleCloseCustomerLists () 
{
	int temp,temp2 ;

	// Shuffling the list of close customers for each customer
	for (int i=nbDepots ; i < nbClients + nbDepots ; i++)
	{
		for (int a1 = 0 ; a1 < (int)cli[i].closeCustomers.size()-1 ; a1++ )
		{
			temp2 = a1 + rand() % ((int)cli[i].closeCustomers.size() - a1) ;
			temp =  cli[i].closeCustomers[a1] ;
			cli[i].closeCustomers[a1] = cli[i].closeCustomers[temp2];
			cli[i].closeCustomers[temp2] = temp ;
		}
	}

	for (int i=nbDepots ; i < nbClients + nbDepots ; i++)
	{
		for (int a1 = 0 ; a1 < (int)cli[i].inverseCloseCustomers.size()-1 ; a1++ )
		{
			temp2 = a1 + rand() % ((int)cli[i].inverseCloseCustomers.size() - a1) ;
			temp =  cli[i].inverseCloseCustomers[a1] ;
			cli[i].inverseCloseCustomers[a1] = cli[i].inverseCloseCustomers[temp2];
			cli[i].inverseCloseCustomers[temp2] = temp ;
		}
	}
}

void Params::ar_computeClientToClientDistances()
{
	// Then, we would still like to include some distance information between services.
	// The distance between two services is the minimum distance between one mode of each service
	// This is used by the granular search
	double myDistanceMin ;
	double myTemp ;

	cluster_to_cluster_cost_lb = new double * [nbClients + nbDepots];
	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		cluster_to_cluster_cost_lb[i] = new double [nbClients + nbDepots];
		for (int j=0 ; j < nbClients + nbDepots ; j++)
		{
			myDistanceMin = 1.e20 ;
			for (int ii=0 ; ii < (int)cli[i].ar_Nodes.size() ; ii++)
			{
				for (int jj=0 ; jj < (int)cli[j].ar_Nodes.size() ; jj++)
				{
					myTemp = cost_mat[cli[i].ar_Nodes[ii]][cli[j].ar_Nodes[jj]] ;
					if (myTemp < myDistanceMin)
						myDistanceMin = myTemp ;
				}
			}
			cluster_to_cluster_cost_lb[i][j] = myDistanceMin ;
		}
	}
}
