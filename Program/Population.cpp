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

#include "Population.h"

Population::Population(Params * params) : params(params)
{
	Individu * randomIndiv ;
	feas_subpop = new SubPop();
	infeas_subpop = new SubPop();
	feas_subpop->nbIndiv = 0 ;
	infeas_subpop->nbIndiv = 0 ;
	double temp, temp2 ;
	bool feasibleFound = false ;

	// Create the trainer
	trainer = new Individu (params,true) ;
	delete trainer->localSearch ;
	trainer->localSearch = new LocalSearch(params,trainer) ; // Initialize the LS structure

	// Creating the initial populations
	for (int i=0 ; i < params->mu && (!params->isSearchingFeasible || !feasibleFound) ; i++ )
	{
		randomIndiv = new Individu (params,true);
		education(randomIndiv);
		addIndividu(randomIndiv) ;
		updateNbFeas(randomIndiv);
		if (!randomIndiv->isValid)
		{
			temp = params->penalityCapa ;
			temp2 = params->penalityLength ;
			params->penalityCapa *= 10 ;
			params->penalityLength *= 10 ;

			trainer->copyIndividual(trainer,randomIndiv);
			trainer->generalSplit();
			trainer->updateLS();
			trainer->localSearch->runSearchTotal();
			trainer->updateIndiv();
			params->penalityCapa = temp ;
			params->penalityLength = temp2 ;
			trainer->generalSplit();
			trainer->copyIndividual(randomIndiv,trainer);
			addIndividu(randomIndiv) ;
		}
		if (randomIndiv->isValid) feasibleFound = true ;
		delete randomIndiv ;
	}

	for (int i=0 ; i < 50 ; i++ )
	{
		if (i%2 == 0) listCapacityFeasible.push_back(true);
		else listCapacityFeasible.push_back(false);
		if (i%2 == 0) listTimeFeasible.push_back(true);
		else listTimeFeasible.push_back(false);
	}

	temp = params->penalityCapa ;
	temp2 = params->penalityLength ;
}

Population::~Population()
{
	int size ;
	if (feas_subpop != NULL)  
	{
		size = (int)feas_subpop->individus.size() ;
		for (int i=0 ; i < size ; i++) delete feas_subpop->individus[i] ;
		delete feas_subpop ;
	}

	if (infeas_subpop != NULL)  
	{
		size = (int)infeas_subpop->individus.size() ;
		for (int i=0 ; i < size ; i++) delete infeas_subpop->individus[i] ;
		delete infeas_subpop ;
	}
	delete trainer ;
}

void Population::evalExtFit(SubPop * pop)
{
	int temp ;
	vector <int> sorted_indices ;
	vector <double> distances ;

	for (int i = 0 ; i < pop->nbIndiv ; i++ )
	{
		sorted_indices.push_back(i) ;
		distances.push_back(pop->individus[i]->avDistanceToClosest(params->nbCountDistMeasure)) ;
	}

	// Ranking the individuals in terms of contribution to diversity
	for (int n = 0 ; n < pop->nbIndiv ; n++ )
	{
		for (int i = 0 ; i < pop->nbIndiv - n - 1 ; i++ )
		{
			if ( distances[sorted_indices[i]] < distances[sorted_indices[i+1]] - 0.000001 )
			{
				temp = sorted_indices[i+1] ;
				sorted_indices[i+1] = sorted_indices[i] ;
				sorted_indices[i] = temp ;
			}
		}
	}

	// Computing the biased fitness
	for (int i = 0 ; i < pop->nbIndiv ; i++ )
	{
		pop->individus[sorted_indices[i]]->divRank = (float)i/(float)(pop->nbIndiv-1) ;
		pop->individus[sorted_indices[i]]->fitRank = (float)sorted_indices[i]/(float)(pop->nbIndiv-1) ;
		pop->individus[sorted_indices[i]]->fitnessEtendu = pop->individus[sorted_indices[i]]->fitRank + ((float)1.0-(float)params->el/(float)pop->nbIndiv) * pop->individus[sorted_indices[i]]->divRank ;
	}
}

int Population::addIndividu (Individu * indiv)
{
	SubPop * subpop ;
	int k, result ;

	if ( indiv->isValid ) subpop = feas_subpop ;
	else subpop = infeas_subpop ;

	result = placeIndividu(subpop,indiv);

	// Keeping only the survivors if the maximum size of the population has been reached
	if (result != -1 && subpop->nbIndiv > params->mu + params->lambda )
	{
		while ( subpop->nbIndiv > params->mu)
		{
			k = selectCompromis(subpop);
			removeIndividu(subpop,k);
		}
	}
	return result ;
}

int Population::addAllIndividus (Population * pop)
{
	Individu * randomIndiv ;
	randomIndiv = new Individu (params,1.0);

	for (int i=0 ; i<pop->feas_subpop->nbIndiv ; i++)
	{
		randomIndiv->copyIndividual(randomIndiv,pop->feas_subpop->individus[i]);
		education(randomIndiv);
		addIndividu(randomIndiv);
	}

	for (int i=0 ; i<pop->infeas_subpop->nbIndiv ; i++)
	{
		randomIndiv->copyIndividual(randomIndiv,pop->infeas_subpop->individus[i]); 
		education(randomIndiv);
		addIndividu(randomIndiv);
	}

	delete randomIndiv ;
	return 1 ;
}

void Population::updateProximity (SubPop * pop, Individu * indiv)
{
	for (int k=0 ; k < pop->nbIndiv ; k++)
	{
		if (pop->individus[k] != indiv) 
		{
			pop->individus[k]->addProche(indiv);
			indiv->addProche(pop->individus[k]);
		}
	}
}

bool Population::fitExist ( SubPop * pop, Individu * indiv )
{
	int count = 0 ;
	double distance = indiv->costSol.evaluation ;
	for (int i=0 ; i < (int)pop->nbIndiv ; i++ )
	{
		if (pop->individus[i]->costSol.evaluation >= (distance - 0.01) && pop->individus[i]->costSol.evaluation <= (distance + 0.01))
			count ++ ;
	}
	if (count <= 1) return false ;
	else return true ;
}

void Population::diversify ()
{
	Individu * randomIndiv ;
	double temp = params->penalityCapa ;
	double temp2 = params->penalityLength ;

	while ( feas_subpop->nbIndiv > (int)(0.3*(double)params->mu))
	{
		delete feas_subpop->individus[feas_subpop->nbIndiv-1] ;
		feas_subpop->individus.pop_back();
		feas_subpop->nbIndiv -- ;
	}

	while ( infeas_subpop->nbIndiv > (int)(0.3*(double)params->mu))
	{
		delete infeas_subpop->individus[infeas_subpop->nbIndiv-1] ;
		infeas_subpop->individus.pop_back();
		infeas_subpop->nbIndiv -- ;
	}

	for (int i=0 ; i < params->mu ; i++ )
	{
		randomIndiv = new Individu (params,true);
		education(randomIndiv);
		addIndividu(randomIndiv) ;
		updateNbFeas(randomIndiv);
		if (!randomIndiv->isValid) 
		{
			temp = params->penalityCapa ;
			temp2 = params->penalityLength ;

			params->penalityCapa *= 50 ;
			params->penalityLength *= 50 ;

			trainer->copyIndividual(trainer,randomIndiv);
			trainer->generalSplit();
			trainer->updateLS();
			trainer->localSearch->runSearchTotal();
			trainer->updateIndiv();
			params->penalityCapa = temp ;
			params->penalityLength = temp2 ;
			trainer->generalSplit();
			trainer->copyIndividual(randomIndiv,trainer);
			addIndividu(randomIndiv) ;
		}
		delete randomIndiv ;
	}
}

int Population::placeIndividu(SubPop * pop, Individu * indiv)
{
	Individu * monIndiv = new Individu (params,false) ;
	monIndiv->copyIndividual (monIndiv , indiv) ;

	bool placed = false ;
	int i = (int)pop->individus.size()-1 ;
	pop->individus.push_back(monIndiv);
	while ( i >= 0 && !placed )
	{
		if (pop->individus[i]->costSol.evaluation >= indiv->costSol.evaluation + 0.001 )
		{
			pop->individus[i+1] = pop->individus[i] ;
			i -- ;
		}
		else
		{
			pop->individus[i+1] = monIndiv ;
			placed = true ;
			pop->nbIndiv ++ ;
			updateProximity (pop, pop->individus[i+1]);
			return i+1 ; // success
		}
	}
	if (!placed)
	{
		pop->individus[0] = monIndiv ;
		placed = true ;
		pop->nbIndiv ++ ;
		updateProximity (pop, pop->individus[0]);
		if (pop == feas_subpop) timeBest = clock();
		return 0 ; // success
	}
	throw string ("erreur placeIndividu") ;
	return -3 ;
}

void Population::removeIndividu(SubPop * pop, int p)
{
	Individu * partant = pop->individus[p];

	// Placing the individual at the end
	for ( int i=p+1 ; i < (int)pop->individus.size() ; i++ )
		pop->individus[i-1] = pop->individus[i] ;

	// Removing it from the population
	pop->individus.pop_back();
	pop->nbIndiv -- ;

	// Removing it from the proximity structures
	for (int i=0 ; i < pop->nbIndiv ; i++ )
		pop->individus[i]->removeProche(partant);

	delete partant ;
}

void Population::validatePen (SubPop * subpop)
{
	Individu * indiv ;

	// Updating Individual Evaluations
	for (int i = 0 ; i < subpop->nbIndiv ; i++)
		subpop->individus[i]->costSol.evaluation = subpop->individus[i]->costSol.distance 
		+ params->penalityCapa * subpop->individus[i]->costSol.capacityViol
		+ params->penalityLength * subpop->individus[i]->costSol.lengthViol ;

	for (int i = 0 ; i < subpop->nbIndiv ; i++)
	{
		for (int j = 0 ; j < subpop->nbIndiv - i - 1; j++)
		{
			if (subpop->individus[j]->costSol.evaluation >= subpop->individus[j+1]->costSol.evaluation + 0.01 )
			{
				indiv = subpop->individus[j] ;
				subpop->individus[j] = subpop->individus[j+1] ;
				subpop->individus[j+1] = indiv ;
			}
		}
	}
}

Individu * Population::getIndividuBinT ()
{
	Individu * individu1 ;
	Individu * individu2 ;
	int place1, place2 ;

	// Picking the first individual in the merge of both subpopulations
	place1 = rand() % (feas_subpop->nbIndiv + infeas_subpop->nbIndiv) ;
	if ( place1 >= feas_subpop->nbIndiv ) 
		individu1 = infeas_subpop->individus[place1 - feas_subpop->nbIndiv] ;
	else 
		individu1 = feas_subpop->individus[place1] ;

	// Picking the second individual in the merge of both subpopulations
	place2 = rand() % (feas_subpop->nbIndiv + infeas_subpop->nbIndiv) ;
	if ( place2 >= feas_subpop->nbIndiv ) 
		individu2 = infeas_subpop->individus[place2 - feas_subpop->nbIndiv] ;
	else 
		individu2 = feas_subpop->individus[place2] ;

	evalExtFit(feas_subpop);
	evalExtFit(infeas_subpop);

	// Keeping the best one
	if (individu1->fitnessEtendu < individu2->fitnessEtendu)
		return individu1 ;
	else
		return individu2 ;
}

Individu * Population::getIndividuPourc (int pourcentage)
{
	int place ;
	// Picking the individual in the 25% best of the valide population, if there are individuals in this set
	if ((feas_subpop->nbIndiv*pourcentage)/100 != 0)
	{
		place = rand() % ((feas_subpop->nbIndiv*pourcentage)/100);
		return feas_subpop->individus[place] ;
	}
	// Picking the individual in the 25% best of the invalide population, if there are individuals in this set
	else if ((infeas_subpop->nbIndiv*pourcentage)/100 != 0)
	{
		place = rand() % ((infeas_subpop->nbIndiv*pourcentage)/100);
		return infeas_subpop->individus[place] ;
	}
	else // If everything fails
	{
		throw string ("ERROR SELECTION POURC") ;
		return NULL ;
	}
}

Individu * Population::getBestFeasibleIndividual ()
{
	if (feas_subpop->nbIndiv != 0) return feas_subpop->individus[0] ;
	else return NULL ;
}

Individu * Population::getBestInfeasibleIndividual ()
{
	if (infeas_subpop->nbIndiv != 0) return infeas_subpop->individus[0] ;
	else return NULL ;
}

bool Population::solutionChecker( vector < vector < int  > > & allRoutes,  vector < vector <int > > & allRoutesNodes, double expectedCost, double expectedMaxRoute)
{
	double totalCost = 0 ;
	double routeCost ;
	double maxRouteLength = 0 ;
	double totalLoad ;
	int cour ;

	// Verification of the load constraints
  for (int r=0 ; r < (int)allRoutes.size() ; r++)
  {
    // For each route
    if (allRoutes[r][0] >= params->nbDepots)
    {
      cout << "SOLUTION CHECKER: No depot at beginning of route" << endl ;
      return false ;
    }

    if (allRoutes[r][allRoutes[r].size()-1] >= params->nbDepots)
    {
      cout << "SOLUTION CHECKER: No depot at end of route" << endl ;
      return false ;
    }

    totalLoad = 0 ;
    for (int i=0 ; i < (int)allRoutes[r].size(); i++ )
    {
      cour = allRoutes[r][i] ;
      totalLoad += params->cli[cour].demand;
    }
    if (totalLoad < 0)
    {
      cout << "SOLUTION CHECKER: Issue with the load" << endl ;
      return false ;
    }

    else if (totalLoad > params->ordreVehicules[r].vehicleCapacity + 0.0001)
    {
      cout << "SOLUTION CHECKER: Violation of load constraint" << endl ;
      return false ;
    }
  }

	// Verification of the solution cost
  for (int r=0 ; r < (int)allRoutes.size() ; r++)
  {
    // For each route, sum the distances
    routeCost = 0 ;
    
    for (int i=0 ; i < (int)allRoutesNodes[r].size()-1; i++ )
    {
      int node_idx1 = allRoutesNodes[r][i];
      int node_idx2 = allRoutesNodes[r][i + 1];
      routeCost += params->cost_mat[node_idx1][node_idx2] ;
    }
    if (routeCost > maxRouteLength) // Updating the maximum route cost
      maxRouteLength = routeCost ;
    totalCost += routeCost ; // Summing the costs
  }

	if (totalCost != expectedCost)
	{
		cout << "SOLUTION CHECKER: Cost is not correct" << endl ;
		return false ;
	}

	return true ; // Success
}

double Population::fractionCapacityFeasible ()
{ 
	int count = 0 ;
	for ( list<bool>::iterator it = listCapacityFeasible.begin(); it != listCapacityFeasible.end(); ++it )
		if (*it == true) count ++ ;

	return double(count)/50. ;
}

double Population::fractionTimeFeasible ()
{ 
	int count = 0 ;
	for ( list<bool>::iterator it = listTimeFeasible.begin(); it != listTimeFeasible.end(); ++it )
		if (*it == true) count ++ ;

	return double(count)/50. ;
}

double Population::getDiversity(SubPop * pop)
{
	double total = 0 ;
	int count = 0 ;
	for ( int i=0 ; i < min(pop->nbIndiv,params->mu) ; i++ )
	{
		for (int j=i+1 ; j < min(pop->nbIndiv,params->mu) ; j++ )
		{
			total += pop->individus[i]->distance(pop->individus[j]);
			count ++ ;
		}
	}
	return total / (double)count ;
} 

double Population::getAverageFeasibleCost ()
{
	double sum_of_penalized_costs = 0 ;
	for (int i=0 ; i < min(feas_subpop->nbIndiv,params->mu) ; i ++)
		sum_of_penalized_costs += feas_subpop->individus[i]->costSol.evaluation ;
	return  sum_of_penalized_costs / min(feas_subpop->nbIndiv,params->mu) ;
}

double Population::getAverageInfeasibleCost ()
{
	double sum_of_penalized_costs = 0 ;
	for (int i=0 ; i <  min(infeas_subpop->nbIndiv,params->mu); i ++)
		sum_of_penalized_costs += infeas_subpop->individus[i]->costSol.evaluation ;
	return  sum_of_penalized_costs / min(infeas_subpop->nbIndiv,params->mu) ;
}

int Population::selectCompromis (SubPop * subpop)
{
	// Selects one individual to be eliminated from the population
	vector <int> sorted_indices ;
	int temp, sortant ;

	evalExtFit(subpop);

	for (int i=0 ; i < subpop->nbIndiv ; i++)
		sorted_indices.push_back(i);

	// Adding a penalty in case of clone (in the objective space or solution space)
	for (int i=1 ; i < subpop->nbIndiv ; i++)
	{
		if (subpop->individus[i]->avDistanceToClosest(1) <= 0.001 ) // in solution space
			subpop->individus[i]->fitnessEtendu += 5 ;
		if (fitExist(subpop,subpop->individus[i])) // in objective space
			subpop->individus[i]->fitnessEtendu += 5 ;	
	}

	// Ranking the elements per extended fitness and selecting out the worst
	for (int n = 0 ; n < subpop->nbIndiv ; n++ )
	{
		for (int i = 0 ; i < subpop->nbIndiv - n - 1 ; i++ )
		{
			if ( subpop->individus[sorted_indices[i]]->fitnessEtendu > subpop->individus[sorted_indices[i+1]]->fitnessEtendu )
			{
				temp = sorted_indices[i+1] ;
				sorted_indices[i+1] = sorted_indices[i] ;
				sorted_indices[i] = temp ;
			}
		}
	}

	sortant = sorted_indices[subpop->nbIndiv-1] ;
	return sortant ;
}

void Population::education(Individu * indiv)
{	
	indiv->copyIndividual(trainer,indiv);
	trainer->generalSplit();
	trainer->updateLS();
	trainer->localSearch->runSearchTotal();
	trainer->updateIndiv();
	indiv->copyIndividual(indiv,trainer);
}

void Population::updateNbFeas (Individu * indiv)
{		
	listCapacityFeasible.push_back(indiv->costSol.capacityViol < 0.001);
	listCapacityFeasible.pop_front();
	listTimeFeasible.push_back(indiv->costSol.lengthViol < 0.001);
	listTimeFeasible.pop_front();
}

void Population::afficheEtat(int nbIter)
{
	// Some traces to observe the status of the population
	cout.precision(8);

	cout << "It " << nbIter << " | Sol " ;

	if (getBestFeasibleIndividual () != NULL)
		cout << getBestFeasibleIndividual()->costSol.distance << " " << getBestFeasibleIndividual()->costSol.routes << " "  ;
	else
		cout << "NO-VALID " ;

	if (getBestInfeasibleIndividual () != NULL)
		cout << getBestInfeasibleIndividual()->costSol.evaluation ;
	else
		cout << "NO-INVALID" ;

	cout << " | Avg " << getAverageFeasibleCost() << " " << getAverageInfeasibleCost()
		<< " | Div " << getDiversity(feas_subpop) << " " << getDiversity(infeas_subpop) << endl
		<< " | Val " << fractionCapacityFeasible() << " " << fractionTimeFeasible()
		<< " | Pen " << params->penalityCapa << " " << params->penalityLength
		<< " | Pop " << feas_subpop->nbIndiv << " " << infeas_subpop->nbIndiv ;
	if (getBestInfeasibleIndividual() != NULL && getBestFeasibleIndividual() == NULL)
		cout << " | Feas : distance " << getBestInfeasibleIndividual()->costSol.distance
		<< " duration " << getBestInfeasibleIndividual()->costSol.lengthViol 
		<< " load " << getBestInfeasibleIndividual()->costSol.capacityViol ;
	cout << endl ;
}
