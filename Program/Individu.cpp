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

#include "Individu.h"

Individu::Individu(Params * params, bool createAllStructures) : params(params)
{
	// Initializing some structures that enable to compute the data on subsequences in the Split algorithm
	myseq = new SeqData(params);
	for (int i=0 ; i<params->nbClients + params->nbDepots +1 ; i++) 
		seq.push_back(new SeqData(params)); 

	vector <int> tempVect ;
	vector <double> tempVectDbl ;
	localSearch = new LocalSearch() ;

	// Initializing the chromosome structures
  for (int j = 0 ; j < params->nombreVehicules ; j++ )
    chromR.push_back(-1);

  successors.resize(params->nbDepots + params->nbClients );
  predecessors.resize(params->nbDepots + params->nbClients );
	for (int i=0 ; i< params->nbDepots + params->nbClients ; i++)
	{
    successors[i] = -1;
    predecessors[i] = -1;
	}

	// If we wish to also create the individual, the local search and Split structures
	if (createAllStructures)
	{
		/* CREATING THE INITIAL INDIVIDUAL */
		vector < vector < int > > tempVect2 ;
		int temp, jj; 
		isFitnessComputed = false ;

		// shuffling
    for (int i = 0 ; i <= (int)chromT.size() - 1 ; i++)
    {
      jj = i + rand() % ((int)chromT.size() - i) ;
      temp = chromT[i] ;
      chromT[i] = chromT[jj] ;
      chromT[jj] = temp ;
    }

		/* CREATING THE SPLIT STRUCTURES */

    for (int i = 0 ; i < params->nombreVehicules + 1; i++)
    {
      pred.push_back(tempVect);
      pred[i].push_back(0);
      for (int j = 0 ; j < (int) params->nbClients + params->nbDepots + 1 ; j++)
        pred[i].push_back(0);
    }

		vector <CostSol> potTemp ;
		CostSol csol ;
		csol.evaluation = 1.e30 ;

    int maxDepots = params->nombreVehicules ;

		for (int i = 0 ; i < maxDepots + 1 ; i++)
		{
			potentiels.push_back(potTemp);
			for (int j = 0 ; j < (int) params->nbClients + params->nbDepots + 1 ; j++)
				potentiels[i].push_back(csol);
		}
		potentiels[0][0].evaluation = 0 ;

		for (int i=0 ;  i< params->nbDepots + params->nbClients ; i++)
		{
			costArcsSplit.push_back(vector<CostSol>());
			for (int j=0 ;  j< params->nbDepots + params->nbClients + 1 ; j++)
				costArcsSplit[i].push_back(CostSol());
		}
	}
}

void Individu::copyIndividual (Individu * destination , Individu * source)
{
	destination->chromT = source->chromT ;
	destination->chromR = source->chromR ;
	destination->costSol.capacityViol = source->costSol.capacityViol ;
	destination->costSol.evaluation = source->costSol.evaluation ;
	destination->costSol.distance = source->costSol.distance ;
	destination->costSol.lengthViol = source->costSol.lengthViol ;
	destination->costSol.routes = source->costSol.routes ;
	destination->isFitnessComputed = source->isFitnessComputed ;
	destination->isFeasible = source->isFeasible ;
	destination->successors = source->successors ;
	destination->predecessors = source->predecessors ;
	destination->nbRoutes = source->nbRoutes ;
	destination->maxRoute = source->maxRoute ;
	destination->potentiels = source->potentiels ;
	destination->pred = source->pred ;
}

Individu::~Individu() 
{ 
	FreeClear(seq);
	delete localSearch ;
	delete myseq ;
}

void Individu::generalSplit()
{
	costSol.evaluation = 0 ;
	costSol.capacityViol = 0 ;
	costSol.distance = 0 ;
	costSol.lengthViol = 0 ;
	costSol.routes = 0 ;

	// performing the Split 
	// we first try the simple split, 
  if (chromT.size() != 0  && splitSimple() == 0) splitLF(); // if its not successful we run the Split with limited fleet

	// Do we have a feasible solution
	if (costSol.capacityViol <= 0.00000001 && costSol.lengthViol <= 0.00000001) 
		isFeasible = true ;
	else 
		isFeasible = false;

	// If split was unsuccessful, we allow more capacity violations
	// Usually, the allowed capacity violation is largely sufficient, and this message would indicate a bug.
	if ( costSol.evaluation > 1.e20 ) 
	{
		cout << "Increasing the capacity violation limit in Split" << endl ;
		params->borne *= 1.1 ;
		generalSplit();
	}

	// A quick test for safety (to avoid any possibility of infinite loops and printouts)
	if (params->borne >= 100.0)
		throw string ("Impossible to Split, most likely a problem of the dataset, aborting the run"); 

	isFitnessComputed = true ;
	measureSol(); // calling a post-processing function which fills all other data structures and verifies the solution cost
	computeSuccessorsAndPredecessors (); // updating the predecessor and successor structure
}

// Simple Split, does not necessarily respect the number of vehicles
int Individu::splitSimple() 
{
	// We will only use the line "0" of the potential and pred data structures
	myseq->initialisation(params->ordreVehicules[0].depotNumber,params,this,false);
	double cost, mydist,mytminex,myloadex;
	int j ;

	// Main code of Split
	for (int i=0 ; i < (int)chromT.size() ; i++ )
	{
		// Compute a route with a single visit
		seq[i]->initialisation(params->ordreVehicules[0].depotNumber,params,this,false);
		j = i ;
		while (j < (int)chromT.size() && seq[j]->load <= params->ordreVehicules[0].vehicleCapacity*params->borne )
		{
			// Extend this route to the next visit
			seq[j+1]->concatOneAfter(seq[j],chromT[j],this);
			cost = seq[0]->evaluation(seq[j+1],myseq,&params->ordreVehicules[0],mydist,mytminex,myloadex); // and evaluate

			// Test if this label is better
			if ( potentiels[0][i].evaluation + cost < potentiels[0][j+1].evaluation )
			{
				potentiels[0][j+1].evaluation = potentiels[0][i].evaluation + cost ;
				potentiels[0][j+1].capacityViol = potentiels[0][i].capacityViol + myloadex ;
				potentiels[0][j+1].distance = potentiels[0][i].distance + mydist ;
				potentiels[0][j+1].lengthViol = potentiels[0][i].lengthViol + mytminex ;
				potentiels[0][j+1].routes = potentiels[0][i].routes + 1 ;
				pred[0][j+1] = i ;
			}
			j++ ;

		}
	}

	// Count the number of routes and see if the soltion is OK
	j = (int)chromT.size() ;
	for (int jj = 0 ; jj < params->nombreVehicules ; jj ++ )
	{
		pred[params->nombreVehicules - jj][j] = pred[0][j] ;
		j = pred[params->nombreVehicules - jj][j] ;
	}

	// If we arrived to the beginning
	if (j == 0)
	{
		// We cumulate the costs
		costSol.capacityViol += potentiels[0][chromT.size()].capacityViol ;
		costSol.evaluation += potentiels[0][chromT.size()].evaluation ;
		costSol.distance += potentiels[0][chromT.size()].distance ;
		costSol.lengthViol += potentiels[0][chromT.size()].lengthViol ;
		costSol.routes += potentiels[0][chromT.size()].routes ;
		initPot(); // we reinitialize the dynamic programming structures for the next use
		return 1 ;
	}
	else
	{
		initPot(); // we reinitialize the dynamic programming structures for the next use
		return 0 ;
	}
}

// Split for problems with limited fleet
// To avoid evaluating several time ("m" times, where m is the number of vehicles) the costs of the arcs 
// (this could be time consuming for some complex VRP variants, for example problems with HOS regulations)
// We already pre-compute and store them in "costArcsSplit"
void Individu::splitLF() 
{ 
	double cost,mydist,mytminex,myloadex;
	myseq->initialisation(params->ordreVehicules[0].depotNumber,params,this,false);

	// preprocessing arc costs
	for (int i=0 ; i < (int)chromT.size() ; i++ )
	{
		seq[i]->initialisation(params->ordreVehicules[0].depotNumber,params,this,false);
		costArcsSplit[i][i].evaluation = seq[i]->evaluation(seq[i],myseq,&params->ordreVehicules[0],mydist,mytminex,myloadex);
		costArcsSplit[i][i].capacityViol =  myloadex ;
		costArcsSplit[i][i].distance = mydist ;
		costArcsSplit[i][i].lengthViol  =  mytminex ;
		for (int j=i ; j < (int)chromT.size() && seq[j]->load <= params->ordreVehicules[0].vehicleCapacity*params->borne ; j++ )
		{
			seq[j+1]->concatOneAfter(seq[j],chromT[j],this);
			costArcsSplit[i][j+1].evaluation = seq[j+1]->evaluation(seq[j+1],myseq,&params->ordreVehicules[0],mydist,mytminex,myloadex);
			costArcsSplit[i][j+1].capacityViol =  myloadex ;
			costArcsSplit[i][j+1].distance = mydist ;
			costArcsSplit[i][j+1].lengthViol  =  mytminex ;
		}
	}


	// for each vehicle
	for (int cam = 0 ; cam < params->nombreVehicules ; cam++)
	{
		int i = 0 ;
		// propagate all labels
		while (i < (int)chromT.size() && potentiels[cam][i].evaluation < 1.e29 )
		{
			cost = costArcsSplit[i][i].evaluation ;
			myloadex = costArcsSplit[i][i].capacityViol ;
			mydist =  costArcsSplit[i][i].distance ;
			mytminex = costArcsSplit[i][i].lengthViol ;
			if ( potentiels[cam][i].evaluation + cost < potentiels[cam+1][i].evaluation )
			{
				potentiels[cam+1][i].evaluation = potentiels[cam][i].evaluation  + cost ;
				potentiels[cam+1][i].capacityViol = potentiels[cam][i].capacityViol + myloadex ;
				potentiels[cam+1][i].distance = potentiels[cam][i].distance + mydist ;
				potentiels[cam+1][i].lengthViol = potentiels[cam][i].lengthViol + mytminex ;
				potentiels[cam+1][i].routes = potentiels[cam][i].routes ;
				pred[cam+1][i] = i ;
			}
			int j = i ;
			while (j < (int)chromT.size() && myloadex <= params->ordreVehicules[cam].vehicleCapacity*(params->borne-1.0) )
			{
				cost = costArcsSplit[i][j+1].evaluation ;
				myloadex = costArcsSplit[i][j+1].capacityViol ;
				mydist =  costArcsSplit[i][j+1].distance ;
				mytminex = costArcsSplit[i][j+1].lengthViol ;
				if ( potentiels[cam][i].evaluation + cost < potentiels[cam+1][j+1].evaluation )
				{
					potentiels[cam+1][j+1].evaluation = potentiels[cam][i].evaluation + cost ;
					potentiels[cam+1][j+1].capacityViol = potentiels[cam][i].capacityViol + myloadex ;
					potentiels[cam+1][j+1].distance = potentiels[cam][i].distance + mydist ;
					potentiels[cam+1][j+1].lengthViol = potentiels[cam][i].lengthViol + mytminex ;
					potentiels[cam+1][j+1].routes = potentiels[cam][i].routes + 1 ;
					pred[cam+1][j+1] = i ;
				}
				j++ ;
			}
			i++ ;
		}
	}

	// we cumulate the costs
	costSol.capacityViol += potentiels[params->nombreVehicules][chromT.size()].capacityViol ;
	costSol.evaluation += potentiels[params->nombreVehicules][chromT.size()].evaluation ;
	costSol.distance += potentiels[params->nombreVehicules][chromT.size()].distance ;
	costSol.lengthViol += potentiels[params->nombreVehicules][chromT.size()].lengthViol ;
	costSol.routes += potentiels[params->nombreVehicules][chromT.size()].routes ;

	// and clean the dynamic programming structures
	initPot();
}

void Individu::measureSol()
{
	int j ;
	nbRoutes = 0 ;
	maxRoute = 0 ;
	double myCost ;
	double totalCost = 0 ;

  // we use the result of Split to fill the other data structures
  j = (int)chromT.size() ;
  for (int jj = 0 ; jj < params->nombreVehicules ; jj ++ )
  {
    if ( j != (int)pred[params->nombreVehicules - jj][j] ) 
    {
      // Beginning and end of this route in the chromosome
      int deb = (int)pred[params->nombreVehicules - jj][j] ;
      int end = j-1 ;

      // Counting the number of routes
      nbRoutes++ ; 
      
      // Constructing the data to evaluate the route
      seq[deb]->initialisation(params->ordreVehicules[0].depotNumber,params,this,false);
      for (int i=deb; i <= end ; i++)
        seq[i+1]->concatOneAfter(seq[i],chromT[i],this);
      myCost = seq[deb]->evaluation(seq[end+1],seq[deb],&params->ordreVehicules[0]);
      
      // Measuring this route to see if its longer than the maxRoute
      if (myCost > maxRoute)
        maxRoute = myCost ;

      // And a quick verification of the solution by summing again the cost (for debugging)
      totalCost += myCost ;
    }
    
    j = (int)pred[params->nombreVehicules - jj][j] ;
    chromR[params->nombreVehicules - jj - 1] = j ;
  }

	// To avoid any problem of numerical imprecision
	costSol.evaluation = costSol.distance + params->penalityCapa * costSol.capacityViol
		+ params->penalityLength * costSol.lengthViol ;
}

void Individu::initPot() 
{
	for (int i = 0 ; i < params->nombreVehicules + 1 ; i++)
		for (size_t j = 0 ; j <= chromT.size() + 1 ; j++)
			potentiels[i][j].evaluation = 1.e30 ;

	potentiels[0][0].evaluation = 0 ;
	potentiels[0][0].capacityViol = 0 ;
	potentiels[0][0].distance = 0 ;
	potentiels[0][0].lengthViol = 0 ;
	potentiels[0][0].routes = 0 ;
}

void Individu::updateLS() 
{
	// Loading the local search structures
	// Warning, Split must have been computed before
	int i,j ;
	ClusterInRoute * myDepot ;
	ClusterInRoute * myDepotFin ;
	ClusterInRoute * myClient ;
	Route * myRoute ;

	// we reinitialize the "shuffledDepotAndCustomerIndices" vector
	localSearch->shuffledDepotAndCustomerIndices.clear() ; 

	// we reset the "isPresent" values to false
	for (i=params->nbDepots ; i < params->nbClients + params->nbDepots ; i++)
		localSearch->clients[i].isPresent = false ;

	// using Split results to load the solution in the LS structure
	j = (int)chromT.size() ;

	for (int jj = 0 ; jj < params->nombreVehicules ; jj ++ )
	{
		i = (int)pred[params->nombreVehicules - jj][j] ;

		myDepot = &localSearch->depots[params->nombreVehicules - jj - 1];
		myDepotFin = &localSearch->depotsFin[params->nombreVehicules - jj - 1];
		myRoute = &localSearch->routes[params->nombreVehicules - jj - 1];

		myDepot->successor = myDepotFin ;
		myDepot->pred = myDepotFin ;
		myDepotFin->successor = myDepot ;
		myDepotFin->pred = myDepot ;

		// a single visit
		if ( j == i+1 )
		{
			myClient = &localSearch->clients[chromT[i]] ;
			myClient->pred = myDepot ;
			myClient->successor = myDepotFin ;
			myClient->route = myRoute ;
			myClient->isPresent = true ;
			localSearch->nodeTestedForEachRoute(myClient->cluster_idx);
			myDepot->successor = myClient ;
			myDepotFin->pred = myClient ;
			localSearch->shuffledDepotAndCustomerIndices.push_back(myClient->cluster_idx);
		}
		else if (j > i+1)
		{
			// at least two visits
			myClient = &localSearch->clients[chromT[i]] ;
			myClient->pred = myDepot ;
			myClient->successor = &localSearch->clients[chromT[i+1]] ;
			myClient->route = myRoute ;
			myClient->isPresent = true ;
			localSearch->nodeTestedForEachRoute(myClient->cluster_idx);
			myDepot->successor = myClient ;
			localSearch->shuffledDepotAndCustomerIndices.push_back(myClient->cluster_idx);

			// information for the end of the route
			myClient = &localSearch->clients[chromT[j-1]] ;
			myClient->pred = &localSearch->clients[chromT[j-2]];
			myClient->successor = myDepotFin ;
			myClient->route = myRoute ;
			myClient->isPresent = true ;
			localSearch->nodeTestedForEachRoute(myClient->cluster_idx);
			myDepotFin->pred = myClient ;
			localSearch->shuffledDepotAndCustomerIndices.push_back(myClient->cluster_idx);

			// and the middle
			for ( int k = (int)i+1 ; k <= j-2 ; k++ )
			{
				myClient = &localSearch->clients[chromT[k]] ;
				myClient->pred = &localSearch->clients[chromT[k-1]] ;
				myClient->successor = &localSearch->clients[chromT[k+1]] ;
				myClient->route =  myRoute ;
				myClient->isPresent = true ;
				localSearch->nodeTestedForEachRoute(myClient->cluster_idx);
				localSearch->shuffledDepotAndCustomerIndices.push_back(myClient->cluster_idx);
			}
		}
		j = i ;
	}
	
	for (i = 0 ; i < params->nombreVehicules ; i++ )
		localSearch->routes[i].updateRouteData(false);
	

	// and we preprocess the route vide structure
  localSearch->setEmptyRoute();
}


void Individu::updateIndiv()
{
	// Now, we go through the LS structure to update the individual (its chromosomes)
	int pos ; 
	vector < Route * > ordreRoutes ;
	ClusterInRoute * node ;

  ordreRoutes.clear();
  for (int ii=0 ; ii < params->nombreVehicules ; ii++)
    ordreRoutes.push_back(&localSearch->routes[ii]);

  pos = 0 ;
  for (int r=0 ; r < (int)ordreRoutes.size() ; r ++)
  {
    node = ordreRoutes[r]->depot->successor ;
    while (!node->isDepot)
    {
      chromT[pos]= node->cluster_idx;
      node = node->successor ;
      pos ++ ;
    }
  }

	// And we apply Split to derive all other structures
	generalSplit();
}

double Individu::distance(Individu * indiv2)
{
	// Hamming distance
	bool isIdentical ;
	double note = 0.0 ;

	for (int j=params->nbDepots ; j < params->nbClients + params->nbDepots ; j++)
	{
		isIdentical = true;

    if (( successors[j] != indiv2->successors[j] || predecessors[j] != indiv2->predecessors[j] )
      && ( predecessors[j] != indiv2->successors[j] || successors[j] != indiv2->predecessors[j] ))
      isIdentical = false ;

		if (!isIdentical)
			note += 1.0 ;
	}

	return ((double)note /(double)(2*params->nbClients)) ;
}

void Individu::computeSuccessorsAndPredecessors ()
{
	int jj ;
	for (int i=0 ; i< params->nbDepots + params->nbClients ; i++)
	{
    successors[i] = -1 ;
    predecessors[i] = -1 ;
	}

  if (chromT.size() != 0)
  {
    for (int i=0 ; i < (int)chromT.size()-1 ; i++)
      successors[chromT[i]] = chromT[i+1];

    for (int i=1 ; i < (int)chromT.size() ; i++)
      predecessors[chromT[i]] = chromT[i-1];

    successors[chromT[chromT.size()-1]] = params->ordreVehicules[0].depotNumber;
    predecessors[chromT[0]] = params->ordreVehicules[0].depotNumber;

    // arranging those which are located at the beginning or end of a route
    for (int i=0 ; i < params->nombreVehicules ; i++)
    {
      jj = chromR[i] ;
      predecessors[chromT[jj]] = params->ordreVehicules[0].depotNumber;
      if (jj != 0)
        successors[chromT[jj-1]] = params->ordreVehicules[0].depotNumber;
    }
  }
}

void Individu::addIndividualProximityInfo(Individu * indiv)
{
	// Adding an individual in the structure of proximity (diversity management procedures)
	list<proxData>::iterator it ;
	proxData data ;
	data.indiv = indiv ;
	data.dist = distance(indiv);

	if (individualProximityInfo.empty()) individualProximityInfo.push_back(data);
	else
	{
		it = individualProximityInfo.begin();
		while ( it != individualProximityInfo.end() && it->dist < data.dist)
			++it ;
		individualProximityInfo.insert(it,data);
	}
}

void Individu::removeIndividualProximityInfo(Individu * indiv)
{
	// Removing an individual in the structure of proximity (diversity management procedures)
	list<proxData>::iterator last = individualProximityInfo.end();
	for (list<proxData>::iterator first = individualProximityInfo.begin() ; first != last ; )
		if (first->indiv == indiv)
			first = individualProximityInfo.erase(first) ;
		else
			++first ;
}

double Individu::avDistanceToClosest(int n) 
{
	// Computing the average distance with the close elements (diversity management)
	double result = 0 ;
	double compte = 0 ;
	list<proxData>::iterator it = individualProximityInfo.begin();

	for (int i=0 ; i<n && it!= individualProximityInfo.end(); i++)
	{
		result += it->dist ;
		compte += 1.0 ;
		++it ;
	}
	return result/compte ;
}
