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

	isILS_general = false ; // Are we running the ILS version of the code
	mu = 25 ; // Population size
	lambda = 40 ; // Number of individuals per generation
	el = 12 ; // Number of elite

	nbCountDistMeasure = 3 ; // Number of close individuals considered in the distance measure (diversity management)
	granularity = 40 ; // Restriction of the LS moves to 40 closest nodes
	minValides = 0.15 ; // Target proportion of feasible solution
	maxValides = 0.20 ; // Target proportion of feasible solution
	penalityCapa = 50 ; // Initial penalties (will evolve during the search)
	penalityLength = 50; // Initial penalties (will evolve during the search)

	// The ELS/ILS requires slightly different parameter setting to get the right number of children and solutions, as specified in Prins 2009
	if (isILS_general) 
	{ mu = 5 ; lambda = 5 ; el = 1 ; minValides = 0.6 ; maxValides = 0.7 ; }
}

void Params::preleveDonnees ()
{
}

bool compPredicate(pairB i, pairB j) 
{ 
	return (i.myparams->timeCost[i.myInt][i.iCour] < i.myparams->timeCost[j.myInt][i.iCour]) ;
}

void Params::calculeStructures () 
{
	// Initializing other structures of the search
	vector < vector < bool > > tempB ;
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
		cli[i].sommetsVoisins.clear();
		cli[i].sommetsVoisinsAvant.clear();
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
			cli[i].sommetsVoisinsAvant.push_back(myVector[j].myInt);
			cli[myVector[j].myInt].sommetsVoisins.push_back(i);
			isCorrelated[myVector[j].myInt][i] = true ;
		}
	}

	for ( int i=0 ; i < nbDepots + nbClients ; i++ )
	{
		cli[i].computeVisitsDyn(nbDays,ancienNbDays);
		cli[i].computeJourSuiv(nbDays,ancienNbDays);
	}
}

void Params::getClient (int i, Client * myCli)
{
	// Reading/Initializing the data for each customer
	string tempstring ;
	pattern p ;
	char * myChars = new char[100]; 
	int t_arc, t_from, t_to, t_qty, t_trav, t_col, t_inv, t_freq ;
	bool isNewService ;

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

	// Setting the pattern information (on which day the customer can be visited)
  p.dep = 0 ;
  p.cost = 0 ;
  p.pat = 1 ;
  myCli->visits.push_back(p);
  myCli->visitsOrigin.push_back(p);

	if (i >= 1) nbTotalServices += myCli->freq ;
	delete [] myChars ;
}

// vc is vehicle capacity. Assume depot cluster is cluster 0 and depot node is node 0
Params::Params(int seedRNG, int nbVeh, double vc, const Ref<const VectorXd> &demand_per_cluster, double time_limit, bool take_first_soln, const std::vector<py::array_t<long>> &cluster_to_node_ptr, const Ref<const MatrixXd> &cost_mat, const Ref<const VectorXl> &node_to_cluster_ptr):type(34), nbVehiculesPerDep(nbVeh), nbDepots(1), isSearchingFeasible(false)
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
	vector <Vehicle> tempI ;
	nbTotalServices = 0 ;
	totalDemand = 0 ;
	multiDepot = false ;
	periodique = false ;

  type = 34;

  int num_clusters = cluster_to_node_ptr.size() - 1; // Not including depot

	// Setting the other instance parameters
	nbDays = 1 ;
	ancienNbDays = nbDays ;
	nbClients = num_clusters;

	// Building the list of vehicles (TODO: why do we have a day zero with zero vehicles and zero capacity?)
	ordreVehicules.push_back(tempI) ;
	nombreVehicules.push_back(0);
	dayCapacity.push_back(0);
	for (int kk=1 ; kk <= nbDays; kk ++)
	{
		ordreVehicules.push_back(tempI) ;
		dayCapacity.push_back(0);
		nombreVehicules.push_back(nbDepots*nbVehiculesPerDep);
		for (int i=0 ; i < nbDepots ; i++)
		{
			for (int j=0 ; j < nbVehiculesPerDep ; j++)
			{
				ordreVehicules[kk].push_back(Vehicle(i,1000000,vc)); // Duration constraint set to a high value
				dayCapacity[kk] += vc ;
			}
		}
	}

	// Reading the list of customers 
	cli = new Client [nbDepots + nbClients] ;
	for (int i = 0 ; i < nbClients + nbDepots ; i++) {
		getClient (i, &cli[i]);
  }

  // Formerly the code from ar_parseOtherLinesNEARP_TP();

  ar_nbArcsDistance = cost_mat.rows();
  ar_Arcs = vector <Arc> (ar_nbArcsDistance) ;
  vector <Arc*> myTempVect = vector <Arc*> (1 + num_clusters, NULL) ;

	/* FIRST SPECIFY THE DEPOT CHARACTERISTICS */
	ar_Arcs[0].indexArc = 0 ;
	ar_Arcs[0].cost = 0 ;
	ar_Arcs[0].nb_Turns = 0 ;
	ar_Arcs[0].nodeBegin = 0 ;
	ar_Arcs[0].nodeEnd = 0 ;
	cli[0].ar_Modes.push_back(&ar_Arcs[0]);
	cli[0].ar_nbModes = 1 ;

	/* THEN READ THE CHARACTERISTICS OF THE NODE SERVICES */
	for (int i=1 ; i <= num_clusters; i++)
	{
    cli[i].ar_nodesExtr0 = i ;
    cli[i].demand = demand_per_cluster[i] ;
    totalDemand += cli[i].demand ;
	}

  for (int i=1 ; i < cost_mat.rows(); i++)
	{
		ar_Arcs[i].indexArc = i ;
		ar_Arcs[i].nodeBegin = i ;
		ar_Arcs[i].nodeEnd = i ;
		ar_Arcs[i].cost = 0. ;
	}

	/* BUILDING THE ARC DISTANCE STRUCTURES */

	ar_distanceArcs.resize(ar_nbArcsDistance);
	for (int i=0 ; i < ar_nbArcsDistance ; i++) {
		ar_distanceArcs[i] = std::vector<double>(ar_nbArcsDistance, 1.e20);
  }

	for (int i=0 ; i < ar_nbArcsDistance ; i++) {
		ar_distanceArcs[i][i] = 0 ; 
  }

  // TODO: this comment was in the original code. Is there something that assumes this elsewhere?
	// distance of any connected pair of edges with a connected at the depot location should have a turn penalty of 0

	for (int i=0 ; i < ar_nbArcsDistance ; i++)
	{
		for (int j=0 ; j < ar_nbArcsDistance ; j++)
		{
      ar_distanceArcs[i][j] = cost_mat(i, j);
		}
	}

	/* SPECIFYING THE MODES FOR THE NODE DELIVERIES */
	for (int i=1 ; i <= num_clusters ; i++)
	{
    auto ptr = cluster_to_node_ptr[i].unchecked();
		cli[i].ar_nbModes = cluster_to_node_ptr[i].shape(0) ;
    for (int local_node_idx = 0; local_node_idx < ptr.shape(0); ++local_node_idx) {
      cli[i].ar_Modes.push_back(&ar_Arcs[ptr(local_node_idx)]);
    }
	}

	/* IDENTIFYING THE MAXIMUM NUMBER OF MODES FOR A DELIVERY */
	ar_maxNbModes = 0 ;
	for (int i=0 ; i <= this->nbClients ; i++) {
		if (cli[i].ar_nbModes > ar_maxNbModes) {
			ar_maxNbModes = cli[i].ar_nbModes ;
    }
  }

  ar_computeDistancesArcs();  // For the NEARP-TP, we rely on the distance between arcs (computed in the line graph)

	// Setting the method parameters
	setMethodParams();

	// Computing the other data structures
	calculeStructures();	
}

Params::~Params(void)
{
	for (int i=0 ; i < nbClients + nbDepots ; i++)
		delete [] timeCost[i] ;
	delete [] timeCost ;
	delete [] cli ;
}

void Params::shuffleProches () 
{
	int temp,temp2 ;

	// Shuffling the list of close customers for each customer
	for (int i=nbDepots ; i < nbClients + nbDepots ; i++)
	{
		for (int a1 = 0 ; a1 < (int)cli[i].sommetsVoisins.size()-1 ; a1++ )
		{
			temp2 = a1 + rand() % ((int)cli[i].sommetsVoisins.size() - a1) ;
			temp =  cli[i].sommetsVoisins[a1] ;
			cli[i].sommetsVoisins[a1] = cli[i].sommetsVoisins[temp2];
			cli[i].sommetsVoisins[temp2] = temp ;
		}
	}

	for (int i=nbDepots ; i < nbClients + nbDepots ; i++)
	{
		for (int a1 = 0 ; a1 < (int)cli[i].sommetsVoisinsAvant.size()-1 ; a1++ )
		{
			temp2 = a1 + rand() % ((int)cli[i].sommetsVoisinsAvant.size() - a1) ;
			temp =  cli[i].sommetsVoisinsAvant[a1] ;
			cli[i].sommetsVoisinsAvant[a1] = cli[i].sommetsVoisinsAvant[temp2];
			cli[i].sommetsVoisinsAvant[temp2] = temp ;
		}
	}
}

void Params::ar_InitializeDistanceNodes()
{
	if (ar_NodesNonRequired+ar_NodesRequired  < 0 || ar_NodesNonRequired+ar_NodesRequired  > 1000000)
		throw string("ERROR WHEN READING : Number of nodes has not been correctly read. A very likely cause is the use of the wrong problem type for a given problem instance");

	// Build the ar_distanceNodes data structures
	vector <double> myTemp = vector <double> (ar_NodesNonRequired+ar_NodesRequired+1) ; 
	for (int i=0 ; i <= ar_NodesNonRequired+ar_NodesRequired ; i++)
		myTemp[i] = 1.e20 ;

	ar_distanceNodes.clear();
	for (int i=0 ; i <= ar_NodesNonRequired+ar_NodesRequired ; i++)
		ar_distanceNodes.push_back(myTemp);
}

void Params::ar_computeDistancesNodes()
{
	for (int ii = 1 ; ii <= ar_NodesNonRequired + ar_NodesRequired ; ii++)
		ar_distanceNodes[ii][ii] = 0 ;

	// simple application of the Floyd Warshall algorithm
	for (int k=1 ; k <= ar_NodesNonRequired + ar_NodesRequired ; k++)
	{
		for (int i=1 ; i <= ar_NodesNonRequired + ar_NodesRequired ; i++)
		{
			for (int j=1 ; j <= ar_NodesNonRequired + ar_NodesRequired ; j++)
			{
				if (ar_distanceNodes[i][k] + ar_distanceNodes[k][j] < ar_distanceNodes[i][j])
					ar_distanceNodes[i][j] = ar_distanceNodes[i][k] + ar_distanceNodes[k][j] ;
			}
		}
	}

	// Then, we would still like to include some distance information between services.
	// The distance between two services is the minimum distance between the closest endpoints of the edge
	// This is used by the granular search
	double d ;
	timeCost = new double * [nbClients + nbDepots+1];
	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		timeCost[i] = new double [nbClients + nbDepots+1];
		for (int j=0 ; j < nbClients + nbDepots ; j++)
		{
			d = min(min(ar_distanceNodes[cli[i].ar_nodesExtr0][cli[j].ar_nodesExtr0],
				ar_distanceNodes[cli[i].ar_nodesExtr0][cli[j].ar_nodesExtr1]),
				min(ar_distanceNodes[cli[i].ar_nodesExtr1][cli[j].ar_nodesExtr0],
				ar_distanceNodes[cli[i].ar_nodesExtr1][cli[j].ar_nodesExtr1]));

			timeCost[i][j] = d ;
		}
	}
}

void Params::ar_computeDistancesArcs()
{
	// computes the distance in the line graph (for the turn penalties)
	// simple application of the Floyd Warshall algorithm
	for (int k=0 ; k < ar_nbArcsDistance ; k++)
	{
		for (int i=0 ; i < ar_nbArcsDistance ; i++)
		{
			for (int j=0 ; j < ar_nbArcsDistance  ; j++)
			{
				if (ar_distanceArcs[i][k] + ar_distanceArcs[k][j] < ar_distanceArcs[i][j])
					ar_distanceArcs[i][j] = ar_distanceArcs[i][k] + ar_distanceArcs[k][j] ;
			}
		}
	}

	// Then, we would still like to include some distance information between services.
	// The distance between two services is the minimum distance between one mode of each service
	// This is used by the granular search
	double myDistanceMin ;
	double myTemp ;

	timeCost = new double * [nbClients + nbDepots];
	for (int i=0 ; i < nbClients + nbDepots ; i++)
	{
		timeCost[i] = new double [nbClients + nbDepots];
		for (int j=0 ; j < nbClients + nbDepots ; j++)
		{
			myDistanceMin = 1.e20 ;
			for (int ii=0 ; ii < (int)cli[i].ar_Modes.size() ; ii++)
			{
				for (int jj=0 ; jj < (int)cli[j].ar_Modes.size() ; jj++)
				{
					myTemp = ar_distanceArcs[cli[i].ar_Modes[ii]->indexArc][cli[j].ar_Modes[jj]->indexArc] ;
					if (myTemp < myDistanceMin)
						myDistanceMin = myTemp ;
				}
			}
			timeCost[i][j] = myDistanceMin ;
		}
	}
}
