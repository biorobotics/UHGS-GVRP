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

#include "Genetic.h"

void Genetic::evolve (int maxIterNonProd, int nbRec)
{
  evolveHGA(maxIterNonProd, nbRec);
}

void Genetic::evolveHGA (int maxIterNonProd, int nbRec)
{
	// Main code of the HGA 
	Individu * parent1 ;
	Individu * parent2 ;
	int place, place2 ;
	place2 = 10000 ;
	nbIterNonProd = 1 ;
	nbIter = 0 ;
	string temp ;
	double fitBeforeRepair ;
	CostSol bestSolFeasibility ;
	clock_t debut = clock() ; 

	if (population->getBestFeasibleIndividual() != NULL) bestSolFeasibility = population->getBestFeasibleIndividual()->costSol ;
	else bestSolFeasibility = population->getBestInfeasibleIndividual()->costSol ;
	for (int i=0 ; i<population->infeas_subpop->nbIndiv ; i++)
		if (population->infeas_subpop->individus[i]->costSol.isBetterFeas(bestSolFeasibility)) bestSolFeasibility = population->infeas_subpop->individus[i]->costSol ;

	offspring->localSearch->nbTotalRISinceBeginning = 0 ;

	cout << "| Start of GA | NbNodes : " << params->nbClients << " | NbVehicles : " << params->nbVehiculesPerDep << " | " << endl ;

	while (nbIterNonProd < maxIterNonProd && (clock() - debut <= ticks) && (!params->isSearchingFeasible || population->getBestFeasibleIndividual() == NULL))
	{
    if (population->getBestFeasibleIndividual() != NULL && population->getBestFeasibleIndividual()->costSol.distance < params->cost_thresh) {
      std::cout << "Returning because we found a solution with cost lower than the provided threshold" << std::endl;
      std::cout << "Cost is " << population->getBestFeasibleIndividual()->costSol.distance << std::endl;
      break;
    }
		// CROSSOVER
		parent1 = population->getIndividuBinT(); // Pick two individuals per binary tournament
		parent2 = population->getIndividuBinT(); // Pick two individuals per binary tournament
		offspring->copyIndividual(offspring,parent1); // Put them in adequate data structures
		parent2_copy->copyIndividual(parent2_copy,parent2); // Put them in adequate data structures

    crossOX(); // OX crossover

		// SPLIT
		offspring->generalSplit();

		// LOCAL SEARCH
		offspring->updateLS();
		offspring->localSearch->runSearchTotal();
		offspring->updateIndiv();
		population->updateNbFeas(offspring);
		place = population->addIndividu(offspring) ;

		// POSSIBLE REPAIR
		if (!offspring->isValid) 
		{
			fitBeforeRepair = offspring->costSol.evaluation ;
			if (rand() % 2 == 0) // 50% chance to do repair on an infeasible individual
			{
				repair();
				if (offspring->costSol.evaluation < fitBeforeRepair - 0.01 || offspring->costSol.evaluation > fitBeforeRepair + 0.01 || offspring->isValid) 
					place2 = population->addIndividu(offspring) ;
				if (offspring->isValid)
					place = place2 ;
				else 
					place = min(place,place2);
			}
		}

		// SOME TRACES
		if ( (offspring->isValid && place == 0) || (offspring->costSol.isBetterFeas(bestSolFeasibility) && population->feas_subpop->nbIndiv == 0))
		{	
			if (traces && population->feas_subpop->nbIndiv > 0) 
				cout << "NEW BEST FEASIBLE " << place << " " << population->getBestFeasibleIndividual()->costSol.evaluation << " distance : " << offspring->costSol.distance << " nbRoutes : " << offspring->costSol.routes << " capaViol : " << offspring->costSol.capacityViol << " lengthViol : " << offspring->costSol.lengthViol << endl << endl ;
			if (traces && population->feas_subpop->nbIndiv == 0 ) 
				cout << "NEW BEST INFEASIBLE "<< place << " " << offspring->costSol.evaluation                            << " distance : " << offspring->costSol.distance << " nbRoutes : " << offspring->costSol.routes << " capaViol : " << offspring->costSol.capacityViol << " lengthViol : " << offspring->costSol.lengthViol << endl << endl ;
			if (offspring->costSol.isBetterFeas(bestSolFeasibility)) 
				bestSolFeasibility = offspring->costSol ;
			nbIterNonProd = 1 ; 
		}
		else nbIterNonProd ++ ;


		// DIVERSIFICATION
		if (nbRec > 0 && nbIterNonProd % (maxIterNonProd/3+1) == maxIterNonProd/3) 
		{	
			if (traces) cout << "Diversification" << endl ;
			population->diversify();
		}

		// PENALTY MANAGEMENT
		if (nbIter % 30 == 0) 
			managePenalites () ;

		// MORE TRACES
		if (traces && nbIter % 500 == 0)
		{
			population->afficheEtat(nbIter);
			cout << " | NbTotalMovesLS : " << offspring->localSearch->nbTotalRISinceBeginning << endl ;
			cout << " | interSwap " << offspring->localSearch->nbInterSwap ;
			cout << " | intraSwap " << offspring->localSearch->nbIntraSwap ;
			cout << " | inter2opt " << offspring->localSearch->nbInter2Opt ;
			cout << " | intra2opt " << offspring->localSearch->nbIntra2Opt ;
			cout << " | " << endl ;
			cout << " | CPU Time : " <<  (double)clock()/(double)CLOCKS_PER_SEC << " seconds " << endl ;
			cout << endl ;
		}
		nbIter ++ ;
	}

	// END OF THE ALGORITHM
	if (traces)
	{
		cout << "Time Elapsed : " << clock() << endl ;
		cout << "Number of Iterations : " << nbIter << endl ;
	}
}

void Genetic::repair ()
{
	double temp, temp2  ;

	temp = params->penalityCapa ;
	temp2 = params->penalityLength ;

	// First Tentative of Repair
	params->penalityCapa *= 10 ;
	params->penalityLength *= 10 ;
	offspring->updateLS();
	offspring->localSearch->runSearchTotal();
	offspring->updateIndiv();

	// If the first tentative failed, second tentative with higher penalty
	if (!offspring->isValid) 
	{
		params->penalityCapa *= 10 ;
		params->penalityLength *= 10 ;
		offspring->generalSplit();
		offspring->updateLS();
		offspring->localSearch->runSearchTotal();
		offspring->updateIndiv();
	}

	params->penalityCapa = temp ;
	params->penalityLength = temp2 ;
	offspring->measureSol();
}

void Genetic::managePenalites ()
{
	bool changeDone = false ;
	double fractionCapacityFeasible = population->fractionCapacityFeasible() ;
	double fractionTimeFeasible = population->fractionTimeFeasible() ;

	// if there are not enough feasible solutions
	if ( fractionCapacityFeasible < params->minFeasibleFraction && params->penalityCapa < 5000)
	{
		params->penalityCapa = (double)((float)params->penalityCapa * 1.2) ;
		changeDone = true ;
	}
	else if ( fractionCapacityFeasible > params->maxFeasibleFraction && params->penalityCapa > 0.01)
	{
		params->penalityCapa =  (double)((float)params->penalityCapa * 0.85) ;
		changeDone = true ;
	}

	// if there are too many feasible solutions
	if ( fractionTimeFeasible < params->minFeasibleFraction && params->penalityLength < 5000)
	{
		params->penalityLength =  (double)((float)params->penalityLength * 1.2) ;
		changeDone = true ;
	}
	else if ( fractionTimeFeasible > params->maxFeasibleFraction && params->penalityLength > 0.01)
	{
		params->penalityLength =  (double)((float)params->penalityLength * 0.85) ;
		changeDone = true ;
	}

	if (changeDone) 
		population->validatePen(population->infeas_subpop);
}

void Genetic::crossOX ()
{
	int temp;

	// We pick the beginning and end of the crossover zone
	int debut = rand() % params->nbClients ;
	int fin = rand() % params->nbClients ;
	while (fin == debut && params->nbClients > 1)
		fin = rand() % params->nbClients ;

	// We initialize a little frequency table to know if each customer was placed or not
	for (int i=params->nbDepots ; i < params->nbClients + params->nbDepots ; i++ )
		freqClient[i] = 1 ;

	int j = debut ;
	// we keep the elements from "debut" to "end"
	while ((j % params->nbClients) != ((fin + 1) % params->nbClients))
	{ 
		freqClient[offspring->chromT[j % params->nbClients]] = 0 ;
		j ++ ;
	}

	// We fill the rest of the elements in the order of the second parent
	for (int i=1 ; i <= params->nbClients ; i++)
	{
		temp = parent2_copy->chromT[(fin + i) % params->nbClients] ;
		if (freqClient[temp] == 1)
		{
			offspring->chromT[j % params->nbClients] = temp ;
			j ++ ;
		}
	}
}

Genetic::Genetic(Params * params,Population * population, clock_t ticks, bool traces) : 
ticks(ticks), traces(traces), population(population), params(params)
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i++ )
		freqClient.push_back(params->cli[i].freq);

	// Creating the Individuals that serve to perform the Local Search and other operations
	offspring = new Individu (params, true) ; 
	parent2_copy = new Individu(params, true) ; 
	offspring->localSearch = new LocalSearch(params,offspring) ;
} 

Genetic::~Genetic(void)
{ 
	delete offspring ;
	delete parent2_copy ;
}

