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

#ifndef INDIVIDU_H
#define INDIVIDU_H

#include <vector>
#include <list>
#include <iostream> 
#include "ClusterInRoute.h"
#include "Params.h"
#include "LocalSearch.h"
using namespace std ;
class LocalSearch ;

// Structure used to keep the cost of a solution
struct CostSol {

	// value of the fitness, including the penalties for possible violations of the constraints
	double evaluation ;

	// value of the distance (not including the penalty values)
	double distance ;

	// amount of capacity violation
	double capacityViol ;

	// amout of distance violation
	double lengthViol ;

	// number of routes
	int routes ;


	// Constructor
	CostSol(){
		evaluation = 0 ;
		distance = 0 ;
		capacityViol = 0 ;
		lengthViol = 0 ;
		routes = 0 ;
	}

	// says if "this" is better in terms of feasibility than "c2"
	bool isBetterFeas(CostSol c2){
		bool isBetter = false ;
		bool isDominated = false ;
		if (capacityViol > c2.capacityViol + 0.0001) isDominated = true ;
		if (lengthViol > c2.lengthViol + 0.0001) isDominated = true ;
		if (capacityViol < c2.capacityViol - 0.0001) isBetter = true ;
		if (lengthViol < c2.lengthViol - 0.0001) isBetter = true ;
		return (isBetter && !isDominated);
	}
};

// preliminary declaration, because proxData depends on Individu
class Individu ;

// data about the proximity of an Individual with regards to the others in the population.
struct proxData {
	// Individual
	Individu * indiv ;

	// Its distance to the others
	double dist ;
};

class Individu
{

public:

	// Access to the parameters of the problem
	Params * params ;

	// fitness of the individual
	double fitnessEtendu ;

	// rank in terms of diversity
	float divRank ;

	// rank in terms of cost
	float fitRank ;

	// evaluation of the solution
	CostSol costSol ;

	// number of routes (if computed)
	int nbRoutes ;

	// length of the maximum route (if computed, the computation is done during split)
	double maxRoute ;

	// The giant tour of the individual
	// chromT[j], jth customer in the giant tour
	vector<int> chromT ;

	// Indices of the beginning of the routes (if already computed via Split)
	// chromR[j] -> route j, gives the index in the chromT of the first customer in this route
	vector<int > chromR ;

	// successors[client] gives the next customer
	// Used to compute the Hamming distance between solutions
	// if the customer does not exist, returns -1
	vector<int> successors ;

	// predecessors[client] gives the previous customer
	// Used to compute the Hamming distance between solutions
	// if the customer does not exist, returns -1
	vector<int> predecessors ;

	// computing the successors and predecessors tables
	void computeSuccessorsAndPredecessors();

	// working table for split (dynamic programming for Split)
	vector < vector < CostSol >> potentiels ;

	// table of predecessors (dynamic programming for Split)
	 vector < vector<int > > pred ;

	// storing the cost of the evaluations of arcs (i,j) in the Split graph
	vector < vector < CostSol> > costArcsSplit ;

	// tells if the fitness has been computed
	bool isFitnessComputed ;

	// tells if the individual is feasible
	bool isValid ;

	// measure of distance from "this" to an individual indiv2
	double distance(Individu * indiv2);

	// Individuals ranked by proximity in the population
	list <proxData> individualProximityInfo ;

	// functions to manage the "individualProximityInfo" structure
	void addIndividualProximityInfo(Individu * indiv) ;
	void removeIndividualProximityInfo(Individu * indiv) ;

	// average distance with the n closest individuals
	double avDistanceToClosest(int n) ;

	// Data structure to perform a LocalSearch. 
	// Only some complete individuals, "offspring" for example in Generic.h
	// possess this structure. The others in Population.h are mainly used as containers
	LocalSearch * localSearch ;

	// Data structure for preprocessing information on sequences during the Split algorithm
	vector<SeqData *> seq ;
	SeqData * myseq ;

	// Split function
	// tries first the simple Split without considering the limit on the number of vehicles
	// if the solution does not respect the number of trips, calls the Split with limited fleet.
	void generalSplit();

	// simple Split function based on a shortest path (does not necessary respect the number of vehicles)
	// returns 1 in case of success, 0 otherwise
	int splitSimple() ;

	// Split function for problems with limited fleet, based on a m-shortest path. 
	// (repeating m times the iteration of Bellman on all nodes)
	void splitLF() ;

	// quick function to fill correctly all data structures, once Split has been performed
	void measureSol() ;

	// initialization of the potentials vector for Split
	void initPot() ;

	// updating the LocalSearch structure with the information of the Individu.
	// Warning, Split must have been computed before 
	void updateLS() ;

	// Updating the individual data structures from the local search information
	void updateIndiv() ;

	// copy of an individual in the other
	// Warning, only copies the chromosomes for storage, not all other structures (potentials and LS)
	void copyIndividual (Individu * destination , Individu * source);

	// constructor of a random individual
	// if the flag "createAllStructures" is set to true, all search structures, including the LS are also initialized
	Individu(Params * params, bool createAllStructures);

	//destructor
	~Individu();

};
#endif
