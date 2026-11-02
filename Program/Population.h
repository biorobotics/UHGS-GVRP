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

#ifndef POPULATION_H
#define POPULATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include "ClusterInRoute.h"
#include "Individu.h"

using namespace std ;

// Structure to manage a sub-population (feasible or infeasible solutions)
struct SubPop
{
	// Individuals
	vector <Individu *> individus ;

	// Number of Individuals
	int nbIndiv ;
};

class Population
{
   private:

   // List to remember which of the 50 previous individuals were feasible in terms of load capacity
   list <bool> listCapacityFeasible ;

   // List to remember which of the 50 previous individuals were feasible in terms of max travel time
   list <bool> listTimeFeasible ;

   // Place an individual in the population, such that the population
   // remains approximately sorted from smallest to largest penalized cost (i.e. the costSol.evaluation value)
   // Returns its position
   int placeIndividu (SubPop * pop, Individu * indiv);

   public:

   // Education procedure (LS)
   void education(Individu * indiv);

   // Access to the parameters of the problem
   Params * params ;

   // clock time when the best individual was found
   clock_t timeBest ;

   // Auxiliary data structure (Individual) with all local search data structures
   // To do the LS on a given individual, we simply copy in this individual and run the LS there.
   Individu * trainer;
	   
   // check if there is already a solution with the same fitness in the population
   bool fitExist ( SubPop * pop, Individu * indiv ) ;

   // compute the biased fitness of the individuals in the population
   void evalExtFit(SubPop * pop);

   // add an individual in the population. If population becomes too large,
   // remove member
   int addIndividu (Individu * indiv);

   // add all individuals from another population
   int addAllIndividus (Population * pop);

   // remove an individual in the population (chosen accoding to the biased fitness)
   void removeIndividu(SubPop * pop, int p);
   
   // subprocedure that chooses an individual to be removed
   int selectCompromis (SubPop * subpop) ; 

   // update the table of distance (Hamming distance) between individuals to know their proximity
   void updateProximity (SubPop * pop, Individu * indiv);

   // Diversification procedure (replace a large part of the population by new random solutions)
   void diversify ();

   // Feasible and Infeasible subpopulations
   SubPop * feas_subpop;
   SubPop * infeas_subpop;

   // Get one individual per binary tournament
   Individu * getIndividuBinT ();

   // Get one individual with uniform probability in a percentage of the best
   Individu * getIndividuPourc (int pourcentage);

   // Get best feasible individual
   Individu * getBestFeasibleIndividual ();

   // Get best infeasible individual
   Individu * getBestInfeasibleIndividual ();

   // when the penalty coefficient change, need to recompute properly the fitness of the individuals in the population
   void validatePen (SubPop * subpop);

   //////////////////////////////////////////////////////////

   // Solution check
   // Verifies the cost and feasibility of the solution
   // Using only the instance data and the shortest path data (not relying on the auxiliary data structures)
   bool solutionChecker( vector < vector < int  > > & allRoutes,  vector < vector <int > > & allRoutesNodes, double expectedCost, double expectedMaxRoute);
   
   /* FUNCTIONS FOR REGULAR PRINTOUTS OF THE STATUS OF THE POPULATION */

   // get the fraction of valid individuals with respect to the load constraint
   double fractionCapacityFeasible () ;

   // get the fraction of valid individuals with respect to the time constraint
   double fractionTimeFeasible () ;

   // get the diversity of the population
   double getDiversity(SubPop * pop);

   // get the average cost of a feasible solution in the population
   double getAverageFeasibleCost ();

   // get the average cost of an infeasible solution in the population
   double getAverageInfeasibleCost ();

   // print a small report of the status of the population
   void afficheEtat(int NbIter);

   // update the count of valid individuals
   void updateNbFeas (Individu * indiv);

   // Constructor
   Population(Params * params);

   // Destructor
   ~Population();
};

#endif
