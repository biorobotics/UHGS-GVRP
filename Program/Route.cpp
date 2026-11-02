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

#include "Route.h"

Route::Route(void){}

Route::Route(int customer_or_depot_idx, Noeud * depot, Vehicle * vehicle, Params * params, Individu * indiv) : params(params), individu(indiv), customer_or_depot_idx(customer_or_depot_idx), depot(depot), vehicle(vehicle) 
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i++ )
	{
		nodeAndRouteTested.push_back(false);
	}
}

Route::~Route(void){}

void Route::reverse ()
{
	// Reversing the order of the nodes in the sequence
	if (!depot->successor->isDepot)
	{
		Noeud * temp ;
		Noeud * myDepot = depot ;
		Noeud * myDepotFin = depot->pred ;
		Noeud * current = depot->successor ;

		while ( !current->isDepot )
		{
			temp = current->successor ;
			current->successor = current->pred ;
			current->pred = temp ;
			current = temp ;
		}

		temp = myDepot->successor ;
		myDepot->successor = myDepotFin->pred ;
		myDepotFin->pred = temp ;
		myDepot->successor->pred = myDepot ;
		myDepotFin->pred->successor = myDepotFin ;
	}
}

void Route::updateRouteData (bool isForPrinting)
{
	bool firstIt ;
	int place = 0 ;
	double nbNodes = 1 ;

	// Computing the auxiliary data on any subsequence (0..i), using forward recursion
	Noeud * noeud = depot ;
	noeud->place = place ;
	noeud->seq0_i->initialisation(noeud->customer_or_depot_idx,params,individu,isForPrinting);
	noeud->seqi_0->initialisation(noeud->customer_or_depot_idx,params,individu,false);

	firstIt = true ;
	while (!noeud->isDepot || firstIt)
	{
		firstIt = false ;
		noeud = noeud->successor ;
		nbNodes ++ ;
		place ++ ;
		noeud->place = place ;
		if (isForPrinting)
			noeud->seq0_i->concatOneAfterWithPathTracking(noeud->pred->seq0_i,noeud->customer_or_depot_idx,individu);
		else
			noeud->seq0_i->concatOneAfter(noeud->pred->seq0_i,noeud->customer_or_depot_idx,individu);
		noeud->seqi_0->concatOneBefore(noeud->pred->seqi_0,noeud->customer_or_depot_idx,individu);
	}

	// Computing the auxiliary data on any subsequence (i..n), using backward recursion
	noeud = depot->pred ;
	noeud->seqi_n->initialisation(noeud->customer_or_depot_idx,params,individu,false);
	noeud->seqn_i->initialisation(noeud->customer_or_depot_idx,params,individu,false);

	firstIt = true ;
	while ( !noeud->isDepot || firstIt )
	{
		firstIt = false ;
		noeud = noeud->pred ;
		noeud->seqi_n->concatOneBefore(noeud->successor->seqi_n,noeud->customer_or_depot_idx,individu);
		noeud->seqn_i->concatOneAfter(noeud->successor->seqn_i,noeud->customer_or_depot_idx,individu);
	}

	// Computing the auxiliary data on any subsequence (i..j), using forward recursion
	// To gain a bit of time, we limit this preprocessing to subsequences such that i..j does not contain more than "sizeSD" elements
	// (More intelligent strategies could be used (e.g., the hierarchical approach of Irnich))
	Noeud * noeudi ;
	Noeud * noeudj ;
	noeudi = depot ;
	for (int i=0 ; i <= depot->pred->place ; i++)
	{
		noeudi->seqi_j[0]->initialisation(noeudi->customer_or_depot_idx,params,individu,false);
		noeudj = noeudi->successor ;
		for (int j=1 ; j <= depot->pred->place - i && j < params->sizeSD ; j++)
		{
			noeudi->seqi_j[j]->concatOneAfter(noeudi->seqi_j[j-1],noeudj->customer_or_depot_idx,individu);
			noeudj = noeudj->successor ;
		}
		noeudi = noeudi->successor ;
	}

	noeudi = depot->pred ;
	for (int i=0 ; i <= depot->pred->place ; i++)
	{
		noeudi->seqj_i[0]->initialisation(noeudi->customer_or_depot_idx,params,individu,false);
		noeudj = noeudi->pred ;
		for (int j=1 ; j <= depot->pred->place - i && j < params->sizeSD ; j++)
		{
			noeudj->seqj_i[j]->concatOneAfter(noeudj->successor->seqj_i[j-1],noeudj->customer_or_depot_idx,individu);
			noeudj = noeudj->pred ;
		}
		noeudi = noeudi->pred ;
	}

	// Checking the route feasibility
	double dist ;
	double violLoad ;
	double violDuration ;
	currentRouteCost = depot->pred->seq0_i->evaluation(depot->pred->seq0_i,depot->seq1,this->vehicle,dist,violDuration,violLoad);
	
	if (violDuration < 0.001 && violLoad < 0.001) 
		isFeasible = true ;
	else 
		isFeasible = false ;
}
