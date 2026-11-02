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

Route::Route(int cluster_idx, ClusterInRoute * depot, Vehicle * vehicle, Params * params, Individu * indiv) : params(params), individu(indiv), cluster_idx(cluster_idx), depot(depot), vehicle(vehicle) 
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
		ClusterInRoute * temp ;
		ClusterInRoute * myDepot = depot ;
		ClusterInRoute * myDepotFin = depot->pred ;
		ClusterInRoute * current = depot->successor ;

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
	ClusterInRoute * cluster = depot ;
	cluster->place = place ;
	cluster->seq0_i->initialisation(cluster->cluster_idx,params,individu,isForPrinting);
	cluster->seqi_0->initialisation(cluster->cluster_idx,params,individu,false);

	firstIt = true ;
	while (!cluster->isDepot || firstIt)
	{
		firstIt = false ;
		cluster = cluster->successor ;
		nbNodes ++ ;
		place ++ ;
		cluster->place = place ;
		if (isForPrinting)
			cluster->seq0_i->concatOneAfterWithPathTracking(cluster->pred->seq0_i,cluster->cluster_idx,individu);
		else
			cluster->seq0_i->concatOneAfter(cluster->pred->seq0_i,cluster->cluster_idx,individu);
		cluster->seqi_0->concatOneBefore(cluster->pred->seqi_0,cluster->cluster_idx,individu);
	}

	// Computing the auxiliary data on any subsequence (i..n), using backward recursion
	cluster = depot->pred ;
	cluster->seqi_n->initialisation(cluster->cluster_idx,params,individu,false);
	cluster->seqn_i->initialisation(cluster->cluster_idx,params,individu,false);

	firstIt = true ;
	while ( !cluster->isDepot || firstIt )
	{
		firstIt = false ;
		cluster = cluster->pred ;
		cluster->seqi_n->concatOneBefore(cluster->successor->seqi_n,cluster->cluster_idx,individu);
		cluster->seqn_i->concatOneAfter(cluster->successor->seqn_i,cluster->cluster_idx,individu);
	}

	// Computing the auxiliary data on any subsequence (i..j), using forward recursion
	// To gain a bit of time, we limit this preprocessing to subsequences such that i..j does not contain more than "sizeSD" elements
	// (More intelligent strategies could be used (e.g., the hierarchical approach of Irnich))
	ClusterInRoute * clusteri ;
	ClusterInRoute * clusterj ;
	clusteri = depot ;
	for (int i=0 ; i <= depot->pred->place ; i++)
	{
		clusteri->seqi_j[0]->initialisation(clusteri->cluster_idx,params,individu,false);
		clusterj = clusteri->successor ;
		for (int j=1 ; j <= depot->pred->place - i && j < params->sizeSD ; j++)
		{
			clusteri->seqi_j[j]->concatOneAfter(clusteri->seqi_j[j-1],clusterj->cluster_idx,individu);
			clusterj = clusterj->successor ;
		}
		clusteri = clusteri->successor ;
	}

	clusteri = depot->pred ;
	for (int i=0 ; i <= depot->pred->place ; i++)
	{
		clusteri->seqj_i[0]->initialisation(clusteri->cluster_idx,params,individu,false);
		clusterj = clusteri->pred ;
		for (int j=1 ; j <= depot->pred->place - i && j < params->sizeSD ; j++)
		{
			clusterj->seqj_i[j]->concatOneAfter(clusterj->successor->seqj_i[j-1],clusterj->cluster_idx,individu);
			clusterj = clusterj->pred ;
		}
		clusteri = clusteri->pred ;
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
