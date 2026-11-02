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

#include "LocalSearch.h"
#include "Individu.h"

void LocalSearch::runSearchTotal ()
{
	// Shuffling the order of move evaluations
	params->shuffleCloseCustomerLists();
	shuffleDepotAndCustomerIndices();
	int nbMoves = 0 ;

	/* RUNNING THE LS */
	// RI -- Route improvement
	updateMoves ();
	reinitAllMoves();
  nbMoves += mutation () ;
	nbTotalRISinceBeginning += nbMoves ;
}

int LocalSearch::mutation ()
{
	// Local Search
	// Based on the classic moves (Relocate, Swap, CROSS, 2-Opt and 2-Opt*)
	rechercheTerminee = false ;
	double costBeforeRemoval, costAfterRemoval;
	bool gainWhenRemoving ;
	int moveEffectue = 0;
	int nbMoves = 0;
	bool routeVideTestee ;
	ClusterInRoute * tempClusterInRoute ;
	int size2 ;

	// We search and apply moves until a local minimum is attained
	while (!rechercheTerminee)
	{
		rechercheTerminee = true ;
		moveEffectue = 0 ;

		// For every node U in random order
		for ( int posU = 0 ; posU < (int)shuffledDepotAndCustomerIndices.size() ; posU ++ )
		{
			posU -= moveEffectue ; // We return on the last node if there has been a move
			nbMoves += moveEffectue ;
			moveEffectue = 0 ;
			routeVideTestee = false ;

			clusterU = &clients[shuffledDepotAndCustomerIndices[posU]];
			routeU = clusterU->route ;
			vehicleU = routeU->vehicle ;
			x = clusterU->successor ;

			// In the CARP, some service removal do not reduce the cost of the route
			// In this case, very few moves can improve the solution (only 2-opt variants), and thus SWAP and RELOCATE variants do not need to be tested
			costBeforeRemoval = clusterU->seq0_i->evaluation(routeU->depot->pred->seq0_i,routeU->vehicle);
			costAfterRemoval = clusterU->seq0_i->evaluation(clusterU->pred->seq0_i,clusterU->successor->seqi_n,routeU->vehicle);
			gainWhenRemoving = (costAfterRemoval < costBeforeRemoval - 0.1) ;

			// For every node V in random order
			size2 = (int)clusterU->moves.size();
			for ( int posV = 0 ; posV < size2 && moveEffectue == 0 ; posV ++ )
			{
				clusterV = &clients[clusterU->moves[posV]] ;
				routeV = clusterV->route ;
				vehicleV = routeV->vehicle ;

				// If we have not yet tested the moves involving the node U and the route of node V
				// (This flag is reset to false as soon as there is a modification in the route)
				if (!clusterV->route->nodeAndRouteTested[clusterU->cluster_idx])
				{
					y = clusterV->successor ;
					if (routeV->cluster_idx != routeU->cluster_idx)
					{
						if (moveEffectue != 1)
						{
							tempClusterInRoute = clusterV ;
							clusterV = clusterV->successor ;
							y = clusterV->successor ;
							// Testing Relocate, Swap, CROSS and I-CROSS (limited to 2 customers) of nodeU and nodeV 
							// Case where they are in different routes
							if (gainWhenRemoving) moveEffectue = interRouteGeneralInsert(); 
							clusterV = tempClusterInRoute ;
							y = clusterV->successor ;
						}

						// 2-Opt*
						if (moveEffectue != 1) 
							moveEffectue = interRoute2Opt ();

						// 2-Opt* (second type, where the routes can be reversed)
						if (moveEffectue != 1) 
							moveEffectue = interRoute2OptInv ();
					}
					else
					{
						tempClusterInRoute = clusterV ;
						clusterV = clusterV->successor ;
						y = clusterV->successor ;

						// Testing Relocate, Swap, CROSS and I-CROSS (limited to 2 customers) of nodeU and nodeV 
						// Case where they are in the same route
						if (moveEffectue != 1 && gainWhenRemoving) 
							moveEffectue = intraRouteGeneralInsertDroite();
						clusterV = tempClusterInRoute ;
						y = clusterV->successor ;
					}
				}
			}

			clusterV = clusterU->successor ;
			routeV = clusterV->route ;
			vehicleV = routeV->vehicle ;
			y = clusterV->successor ;
			if (!clusterV->route->nodeAndRouteTested[clusterU->cluster_idx])
			{
				while (moveEffectue != 1 && !clusterV->isDepot)
				{ 
					// Testing 2-Opt between U and V (if the restriction of the granular search allows) 
					if (params->isCorrelated[clusterU->pred->cluster_idx][clusterV->cluster_idx] || params->isCorrelated[clusterU->cluster_idx][clusterV->successor->cluster_idx]) 
						moveEffectue = intraRoute2Opt ();
					clusterV = clusterV->successor ;
					y = clusterV->successor ;
				}
			}

			// Special cases : testing the insertions behind the depot, and the empty routes
			for (int route = 0 ; route < params->nbVehiculesPerDep && moveEffectue == 0 ; route ++)
			{
				clusterV = &depots[route] ;
				routeV = clusterV->route ;
				y = clusterV->successor ;
				if ( (!clusterV->route->nodeAndRouteTested[clusterU->cluster_idx]) && (!y->isDepot || !routeVideTestee))
				{
					if (y->isDepot) routeVideTestee = true ;
					if (routeV != routeU)
					{
						tempClusterInRoute = clusterV ;
						clusterV = depots[route].successor ;
						y = clusterV->successor ;

						// Insertion after the depot, in a different route
						if (gainWhenRemoving && (params->isCorrelated[clusterU->cluster_idx][clusterV->cluster_idx] || params->isCorrelated[clusterU->cluster_idx][y->cluster_idx]) && moveEffectue != 1 ) 
							moveEffectue = interRouteGeneralInsert();

						clusterV = clusterV->route->depot ;
						y = clusterV->successor ;

						// 2-Opt* after the depot
						if (params->isCorrelated[clusterU->pred->cluster_idx][clusterV->cluster_idx] && moveEffectue != 1) 
							moveEffectue = interRoute2Opt ();

						// 2-Opt* after the depot
						if ((params->isCorrelated[x->cluster_idx][y->cluster_idx] || params->isCorrelated[y->cluster_idx][x->cluster_idx]) && moveEffectue != 1) 
							moveEffectue = interRoute2OptInv ();

						clusterV = tempClusterInRoute ;
						y = clusterV->successor ;
					}
					else
					{
						tempClusterInRoute = clusterV ;
						clusterV = depots[route].successor ;
						y = clusterV->successor ;

						// Insertion after the depot, in the same route
						if ((params->isCorrelated[clusterU->cluster_idx][clusterV->cluster_idx] || params->isCorrelated[clusterU->cluster_idx][y->cluster_idx]) && moveEffectue != 1) 
							moveEffectue = intraRouteGeneralInsertDroite();

						clusterV = tempClusterInRoute ;
						y = clusterV->successor ;
					}
				}
			}

			// Say that we have tested the node U with all routes
			if (moveEffectue == 0)
				nodeTestedForEachRoute(clusterU->cluster_idx);
		}
	}
	return nbMoves ;
}

int LocalSearch::interRouteGeneralInsert()
{
	// For a pair of nodes U, V, tests together the Relocate, Swap, and variants of CROSS and I-CROSS limited to two consecutive nodes.
	// Some route evaluations can be gained (about 40%) by doing these computations is a combined manner
	int ibest, jbest ;
	double moveMin ;
	double temp ;
	SeqData * seq = clusterU->seq0_i ;
	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ;

	// This table will keep the result of the moves
	// 0 -> send nothing
	// 1 -> send U
	// 2 -> send U,Unext
	// 3 -> send Unext,U
	for (int i=0 ; i<4 ; i++ )
		for (int j=0 ; j<4 ; j++)
			resultMoves[i][j] = 1.e29 ; 

	// This keeps the current cost
	resultMoves[0][0] = costZero ;

	/* EVALUATION OF MOVE LOWER BOUNDS */
	// We start to compute lower bounds on move evaluations, so that we can discard any moves which are not promising.

	// Can send something if U is not a depot
	if (!clusterU->isDepot)
	{
		// We evaluate the result of the move as a combination of existing subsequences
		resultMoves[1][0] = seq->evaluationLB(clusterU->pred->seq0_i,x->seqi_n,routeU->vehicle);

		// Can we receive something (if V is not a depot)
		// In the following, all the conditionals are set up to avoid moving depots
		if (!clusterV->isDepot)
			resultMoves[1][1] = seq->evaluationLB(clusterU->pred->seq0_i,clusterV->seq1,x->seqi_n,routeU->vehicle);

		// And so on...
		if (!x->isDepot)
		{
			resultMoves[2][0] = seq->evaluationLB(clusterU->pred->seq0_i,x->successor->seqi_n,routeU->vehicle);
			resultMoves[3][0] = resultMoves[2][0] ;
			if (!clusterV->isDepot)
			{
				resultMoves[2][1] = seq->evaluationLB(clusterU->pred->seq0_i,clusterV->seq1,x->successor->seqi_n,routeU->vehicle);
				resultMoves[3][1] = resultMoves[2][1] ;
				if (!y->isDepot)
				{
					resultMoves[2][2] = seq->evaluationLB(clusterU->pred->seq0_i,clusterV->seq12,x->successor->seqi_n,routeU->vehicle);
					resultMoves[2][3] = seq->evaluationLB(clusterU->pred->seq0_i,clusterV->seq21,x->successor->seqi_n,routeU->vehicle);
					resultMoves[3][2] = resultMoves[2][2];
					resultMoves[3][3] = resultMoves[2][3];
				}
			}
		}
	}

	if (!clusterU->isDepot)
	{
		resultMoves[1][0] += seq->evaluationLB(clusterV->pred->seq0_i,clusterU->seq1,clusterV->seqi_n,routeV->vehicle);
		if (!x->isDepot)
		{
			resultMoves[2][0] += seq->evaluationLB(clusterV->pred->seq0_i,clusterU->seq12,clusterV->seqi_n,routeV->vehicle);
			resultMoves[3][0] += seq->evaluationLB(clusterV->pred->seq0_i,clusterU->seq21,clusterV->seqi_n,routeV->vehicle);
		}
	}

	if (!clusterV->isDepot)
	{
		if (!clusterU->isDepot)
		{
			resultMoves[1][1] += seq->evaluationLB(clusterV->pred->seq0_i,clusterU->seq1,y->seqi_n,routeV->vehicle);
			if (!x->isDepot)
			{
				resultMoves[2][1] += seq->evaluationLB(clusterV->pred->seq0_i,clusterU->seq12,y->seqi_n,routeV->vehicle);
				resultMoves[3][1] += seq->evaluationLB(clusterV->pred->seq0_i,clusterU->seq21,y->seqi_n,routeV->vehicle);
			}
		}
		if (!y->isDepot && !clusterU->isDepot && !x->isDepot)
		{
			temp = seq->evaluationLB(clusterV->pred->seq0_i,clusterU->seq12,y->successor->seqi_n,routeV->vehicle);
			resultMoves[2][2] += temp ;
			resultMoves[2][3] += temp ;
			temp = seq->evaluationLB(clusterV->pred->seq0_i,clusterU->seq21,y->successor->seqi_n,routeV->vehicle);
			resultMoves[3][2] += temp ;
			resultMoves[3][3] += temp ;
		}
	}

	/* WE IDENTIFY WHICH MOVES CAN BE FILTERED OUT USING THE LOWER BOUND */

	for (int i=0 ; i<4; i++)
	{
		for (int j=0 ; j < 4 ; j++)
		{
			shouldBeTested[i][j] = (resultMoves[i][j] < costZero) ;
		}
	}
	shouldBeTested[0][0] = true ;
	resultMoves[0][0] = costZero ;

	/* AND NOW WE TEST THE MOVES THAT HAVE A CHANCE TO BE IMPROVING */
	// Exactly the same code as previously, but using "seq->evaluation" instead of "seq->evaluationLB"

	if (!clusterU->isDepot)
	{
		if (shouldBeTested[1][0]) resultMoves[1][0] = seq->evaluation(clusterU->pred->seq0_i,x->seqi_n,routeU->vehicle);

		if (!clusterV->isDepot && shouldBeTested[1][1])
			resultMoves[1][1] = seq->evaluation(clusterU->pred->seq0_i,clusterV->seq1,x->seqi_n,routeU->vehicle);

		if (!x->isDepot)
		{
			if (shouldBeTested[2][0] || shouldBeTested[3][0]) 
			{
				resultMoves[2][0] = seq->evaluation(clusterU->pred->seq0_i,x->successor->seqi_n,routeU->vehicle);
				resultMoves[3][0] = resultMoves[2][0] ;
			}
			if (!clusterV->isDepot)
			{
				if (shouldBeTested[2][1] || shouldBeTested[3][1])
				{
					resultMoves[2][1] = seq->evaluation(clusterU->pred->seq0_i,clusterV->seq1,x->successor->seqi_n,routeU->vehicle);
					resultMoves[3][1] = resultMoves[2][1] ;
				}
				if (!y->isDepot)
				{
					if (shouldBeTested[2][2] || shouldBeTested[3][2])
					{
						resultMoves[2][2] = seq->evaluation(clusterU->pred->seq0_i,clusterV->seq12,x->successor->seqi_n,routeU->vehicle);
						resultMoves[3][2] = resultMoves[2][2];
					}
					if (shouldBeTested[2][3] || shouldBeTested[3][3])
					{
						resultMoves[2][3] = seq->evaluation(clusterU->pred->seq0_i,clusterV->seq21,x->successor->seqi_n,routeU->vehicle);
						resultMoves[3][3] = resultMoves[2][3];
					}
				}
			}
		}
	}

	if (!clusterU->isDepot)
	{
		if (shouldBeTested[1][0]) resultMoves[1][0] += seq->evaluation(clusterV->pred->seq0_i,clusterU->seq1,clusterV->seqi_n,routeV->vehicle);
		if (!x->isDepot)
		{
			if (shouldBeTested[2][0]) resultMoves[2][0] += seq->evaluation(clusterV->pred->seq0_i,clusterU->seq12,clusterV->seqi_n,routeV->vehicle);
			if (shouldBeTested[3][0]) resultMoves[3][0] += seq->evaluation(clusterV->pred->seq0_i,clusterU->seq21,clusterV->seqi_n,routeV->vehicle);
		}
	}

	if (!clusterV->isDepot)
	{
		if (!clusterU->isDepot)
		{
			if (shouldBeTested[1][1]) resultMoves[1][1] += seq->evaluation(clusterV->pred->seq0_i,clusterU->seq1,y->seqi_n,routeV->vehicle);
			if (!x->isDepot)
			{
				if (shouldBeTested[2][1]) resultMoves[2][1] += seq->evaluation(clusterV->pred->seq0_i,clusterU->seq12,y->seqi_n,routeV->vehicle);
				if (shouldBeTested[3][1]) resultMoves[3][1] += seq->evaluation(clusterV->pred->seq0_i,clusterU->seq21,y->seqi_n,routeV->vehicle);
			}
		}
		if (!y->isDepot && !clusterU->isDepot && !x->isDepot)
		{
			if (shouldBeTested[2][2] || shouldBeTested[2][3])
			{
				temp = seq->evaluation(clusterV->pred->seq0_i,clusterU->seq12,y->successor->seqi_n,routeV->vehicle);
				resultMoves[2][2] += temp ;
				resultMoves[2][3] += temp ;
			}
			if (shouldBeTested[3][2] || shouldBeTested[3][3])
			{
				temp = seq->evaluation(clusterV->pred->seq0_i,clusterU->seq21,y->successor->seqi_n,routeV->vehicle);
				resultMoves[3][2] += temp ;
				resultMoves[3][3] += temp ;
			}
		}
	}

	// We identify the best move among relocate, swap, CROSS and I-CROSS between the node pair U,V
	ibest = 0 ; jbest = 0 ;
	moveMin = 1.e30 ;
	for (int i=0 ; i<4; i++ )
	{
		for (int j=0 ; j < 4 ; j++)
		{
			if (shouldBeTested[i][j] && resultMoves[i][j] < moveMin - EPSILON_LS )
			{
				moveMin = resultMoves[i][j] ;
				ibest = i ;
				jbest = j ;
			}
		}
	}

	// If no improving move between U and V, we return
	if ( ibest == 0 && jbest == 0) 
		return 0 ;

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	ClusterInRoute * placeV = clusterV->pred ;
	ClusterInRoute * placeU = clusterU->pred ;
	reinitMoves(placeU->route); // Say that all moves involving this route must be tested again
	reinitMoves(placeV->route); // Say that all moves involving this route must be tested again

	// Moving the nodes in the solution structure
	if ( ibest == 1 || ibest == 2) insertClusterInRoute(clusterU,placeV);
	if ( ibest == 2) insertClusterInRoute(x,clusterU);
	if ( ibest == 3 ) { insertClusterInRoute(x,placeV); insertClusterInRoute(clusterU,x); }

	if ( jbest == 1 || jbest == 2) insertClusterInRoute(clusterV,placeU);
	if ( jbest == 2) insertClusterInRoute(y,clusterV);
	if ( jbest == 3 ) { insertClusterInRoute(y,placeU); insertClusterInRoute(clusterV,y); }

	// Update the pre-processed data on the subsequences of the route
	placeU->route->updateRouteData(false);
	placeV->route->updateRouteData(false);
	setRouteVide(); // Keep a pointer on the first empty route

	rechercheTerminee = false ; // Not finished the search
	nbInterSwap ++ ;
	return 1 ; // Return Success
}

int LocalSearch::interRoute2Opt ()
{
	// Testing 2-Opt* between U and V
	double cost ;
	SeqData * seq = clusterU->seq0_i ;

	if  (routeU->depot->cluster_idx != routeV->depot->cluster_idx || (clusterU->pred->isDepot && clusterV->isDepot)) 
		return 0 ; // Cannot do a 2-Opt* if the depot is placed in the wrong way

	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ; // Cost before the move

	// Lower bound on the move value
	cost = seq->evaluationLB(clusterU->pred->seq0_i,y->seqi_n,routeU->vehicle) + seq->evaluationLB(clusterV->seq0_i,clusterU->seqi_n,routeV->vehicle) - costZero ;
	if ( cost  > -EPSILON_LS ) 
		return 0 ; // Exit if no chance of improvement

	// Exact move evaluation
	cost = seq->evaluation(clusterU->pred->seq0_i,y->seqi_n,routeU->vehicle) + seq->evaluation(clusterV->seq0_i,clusterU->seqi_n,routeV->vehicle) - costZero ;
	if ( cost  > -EPSILON_LS ) 
		return 0 ; 

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	reinitMoves(routeU);  // Say that all moves involving this route must be tested again
	reinitMoves(routeV);  // Say that all moves involving this route must be tested again
	ClusterInRoute * tempU = clusterU ;
	clusterU = clusterU->pred ;
	x = clusterU->successor ;

	// Updating the solution
	ClusterInRoute * count ;
	ClusterInRoute * depotU = routeU->depot ;
	ClusterInRoute * depotV = routeV->depot ;
	ClusterInRoute * depotUFin = depotU->pred ;
	ClusterInRoute * depotVFin = depotV->pred ;
	ClusterInRoute * depotUpred = depotUFin->pred ;

	count = y ;
	while ( !count->isDepot )
	{
		count->route = routeU ;
		count = count->successor ;
	}

	count = x ;
	while ( !count->isDepot )
	{
		count->route = routeV ;
		count = count->successor ;
	}

	clusterU->successor = y ;
	y->pred = clusterU ;
	clusterV->successor = x ;
	x->pred = clusterV ;

	if (x->isDepot)
	{
		depotUFin->pred = depotVFin->pred ;
		depotUFin->pred->successor = depotUFin ;
		clusterV->successor = depotVFin ;
		depotVFin->pred = clusterV ;
	}
	else
	{
		depotUFin->pred = depotVFin->pred ;
		depotUFin->pred->successor = depotUFin ;
		depotVFin->pred = depotUpred ;
		depotVFin->pred->successor = depotVFin ;
	}

	// Update the pre-processed data on the subsequences of the route
	routeU->updateRouteData(false);
	routeV->updateRouteData(false);
	setRouteVide(); // Keep a pointer on the first empty route

	rechercheTerminee = false ; // Not finished the search
	nbInter2Opt ++ ;
	clusterU = tempU ;
	x = clusterU->successor ;
	return 1 ; // Return Success
}


int LocalSearch::interRoute2OptInv()
{
	// 2-Opt* with route inversions
	SeqData * seq = clusterU->seq0_i ;
	double cost ;
	double costTemp ;
	double costTempReverse = 1.e30 ;
	bool reverseRouteU, reverseRouteV ; 
	if (clusterU->route == clusterV->route) { return 0 ; }
	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ;

	// Compute a lower bound on the value of the move
	cost = - costZero ;
	costTemp = seq->evaluationLB(y->seqn_i,x->seqi_n,routeU->vehicle) ;
	costTempReverse = seq->evaluationLB(x->seqn_i,y->seqi_n,routeU->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteU = false ; }
	else  { cost += costTempReverse ; reverseRouteU = true ; }

	costTemp = seq->evaluationLB(clusterV->seq0_i,clusterU->seqi_0,routeV->vehicle) ;
	costTempReverse = seq->evaluationLB(clusterU->seq0_i,clusterV->seqi_0,routeV->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteV = false ; }
	else  { cost += costTempReverse ; reverseRouteV = true ; }

	if ( cost  > -EPSILON_LS ) 
		return 0 ;  // Exit if no chance of improvement

	// Test the real move cost
	cost = - costZero ;
	costTemp = seq->evaluation(y->seqn_i,x->seqi_n,routeU->vehicle) ;
	costTempReverse = seq->evaluation(x->seqn_i,y->seqi_n,routeU->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteU = false ; }
	else  { cost += costTempReverse ; reverseRouteU = true ; }

	costTemp = seq->evaluation(clusterV->seq0_i,clusterU->seqi_0,routeV->vehicle) ;
	costTempReverse = seq->evaluation(clusterU->seq0_i,clusterV->seqi_0,routeV->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteV = false ; }
	else  { cost += costTempReverse ; reverseRouteV = true ; }

	if ( cost  > -EPSILON_LS ) 
		return 0 ;

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	reinitMoves(clusterU->route);  // Say that all moves involving this route must be tested again
	reinitMoves(clusterV->route);  // Say that all moves involving this route must be tested again
	ClusterInRoute * depotU = routeU->depot ;
	ClusterInRoute * depotV = routeV->depot ;
	ClusterInRoute * depotUFin = routeU->depot->pred ;
	ClusterInRoute * depotVFin = routeV->depot->pred ;
	ClusterInRoute * depotUSuiv = depotU->successor ;

	// Update the solution
	ClusterInRoute * temp ;
	ClusterInRoute * yy = y ;
	ClusterInRoute * uu = clusterU ;

	while ( !yy->isDepot )
	{
		temp = yy->successor ;
		yy->successor = yy->pred ;
		yy->pred = temp ;
		yy->route = routeU ;
		yy = temp ;
	}

	while ( !uu->isDepot )
	{
		temp = uu->pred ;
		uu->pred = uu->successor ;
		uu->successor = temp ;
		uu->route = routeV ;
		uu = temp ;
	}

	clusterV->successor = clusterU ;
	clusterU->pred = clusterV ;
	y->successor = x ;
	x->pred = y ;

	if (y->isDepot)
	{
		depotVFin->successor = depotV ;
		depotVFin->pred = depotUSuiv ;
		depotVFin->pred->successor = depotVFin ;
		depotU->successor = x ;
		x->pred = depotU ;
	}
	else if (clusterU->isDepot )
	{
		depotU->successor = depotVFin->pred ;
		depotU->successor->pred = depotU ;
		depotU->pred = depotUFin ;
		depotVFin->pred = clusterV ;
		clusterV->successor = depotVFin ;
	}
	else
	{
		depotU->successor = depotVFin->pred ;
		depotU->successor->pred = depotU ;
		depotVFin->pred = depotUSuiv ;
		depotVFin->pred->successor = depotVFin ;
	}

	// Reverse if needed
	if (reverseRouteU) routeU->reverse();
	if (reverseRouteV) routeV->reverse();

	// Update the pre-processed data on the subsequences of the route
	routeU->updateRouteData(false);
	routeV->updateRouteData(false);
	setRouteVide(); // Keep a pointer on the first empty route

	rechercheTerminee = false ; // Not finished the search
	nbInter2Opt ++ ;
	return 1 ; // Return Success
}

int LocalSearch::intraRouteGeneralInsertDroite ()
{
	// For a pair of nodes U, V, IN THE SAME ROUTE, tests together the Relocate, Swap, and variants of CROSS and I-CROSS limited to two consecutive nodes.
	ClusterInRoute * tempU = clusterU ;
	ClusterInRoute * tempV = clusterV ;
	bool turned = false ;

	int decalage = clusterV->place - clusterU->place ; // Decalage is the difference of index between U and V
	if (decalage >= -1 && decalage <= 1) return 0 ; // This means that U and V are already consecutive, testing these moves is useless

	// If decalage < 0, this means that V is on the left of U in the route
	// We don't want to deal with this case, so we simply inverse the notations of U and V 
	// And we deal in the following with only the case where V is on the right
	if ( decalage < 0 )
	{
		clusterU = tempV ;
		clusterV = tempU ;
		x = clusterU->successor ;
		y = clusterV->successor ;
		decalage = -decalage ;
		turned = true ;
	}

	double moveMin ;
	int ibest, jbest ;
	SeqData * seq = clusterU->seq0_i ;
	ClusterInRoute * Upred = clusterU->pred ;
	ClusterInRoute * placeU = clusterU->pred ;
	ClusterInRoute * placeV = clusterV->pred ;
	ClusterInRoute * xsuccessor = x->successor ;
	ClusterInRoute * ysuccessor = y->successor ;
	double costZero = routeU->currentRouteCost ;

	// This table will keep the result of the moves
	// 0 -> send nothing
	// 1 -> send U
	// 2 -> send U,Unext
	// 3 -> send Unext,U
	for (int i=0 ; i<4 ; i++ )
		for (int j=0 ; j<4 ; j++ )
			resultMoves[i][j] = 1.e29 ;

	resultMoves[0][0] = costZero ;

	if (decalage >= 3) // General case
	{
		if (!clusterV->isDepot && turned) // relocate of V before U (V cannot be a depot, but U can be)
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			myseqs.push_back(clusterV->seq1);
			// The intra-route moves are the only moves which can involve general subsequences (i,j) in the middle of the route
			// To reduce a bit the pre-processing, we keep only sequences of limited size (<= 15), even if this involves to concatenate more pieces
			// The concatenation of the good number of pieces is done by "addSeqDataInPieces"
			addSeqDataInPieces(clusterU,clusterV->place-1-clusterU->place);
			myseqs.push_back(y->seqi_n);
			resultMoves[0][1] = seq->evaluationLB(myseqs,routeU->vehicle);
			if (!y->isDepot) // exchange U and (V,Y or Y,V)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(clusterV->seq12);
				addSeqDataInPieces(clusterU,clusterV->place-1-clusterU->place);
				myseqs.push_back(ysuccessor->seqi_n);
				resultMoves[0][2] = seq->evaluationLB(myseqs,routeU->vehicle);
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(clusterV->seq21);
				addSeqDataInPieces(clusterU,clusterV->place-1-clusterU->place);
				myseqs.push_back(ysuccessor->seqi_n);
				resultMoves[0][3] = seq->evaluationLB(myseqs,routeU->vehicle);
			}
		}

		if (!turned)
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			addSeqDataInPieces(x,clusterV->place-1-x->place);
			myseqs.push_back(clusterU->seq1);
			myseqs.push_back(clusterV->seqi_n);
			resultMoves[1][0] = seq->evaluationLB(myseqs,routeU->vehicle);
		}
		if (!clusterV->isDepot) // exchange U and V
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			myseqs.push_back(clusterV->seq1);
			addSeqDataInPieces(x,clusterV->place-1-x->place);
			myseqs.push_back(clusterU->seq1);
			myseqs.push_back(y->seqi_n);
			resultMoves[1][1] = seq->evaluationLB(myseqs,routeU->vehicle);
			if (!y->isDepot && turned) // exchange U and (V,Y or Y,V)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(clusterV->seq12);
				addSeqDataInPieces(x,clusterV->place-1-x->place);
				myseqs.push_back(clusterU->seq1);
				myseqs.push_back(ysuccessor->seqi_n);
				resultMoves[1][2] = seq->evaluationLB(myseqs,routeU->vehicle);

				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(clusterV->seq21);
				addSeqDataInPieces(x,clusterV->place-1-x->place);
				myseqs.push_back(clusterU->seq1);
				myseqs.push_back(ysuccessor->seqi_n);
				resultMoves[1][3] = seq->evaluationLB(myseqs,routeU->vehicle);
			}
		}

		if (!x->isDepot)
		{
			if (!turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
				myseqs.push_back(clusterU->seq12);
				myseqs.push_back(clusterV->seqi_n);
				resultMoves[2][0] = seq->evaluationLB(myseqs,routeU->vehicle);
			}

			if (!clusterV->isDepot) // exchange U,X and V
			{
				if (!turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq1);
					addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
					myseqs.push_back(clusterU->seq12);
					myseqs.push_back(y->seqi_n);
					resultMoves[2][1] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
				if (!y->isDepot) // exchange U,X and (V,Y or Y,V)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq12);
					addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
					myseqs.push_back(clusterU->seq12);
					myseqs.push_back(ysuccessor->seqi_n);
					resultMoves[2][2] = seq->evaluationLB(myseqs,routeU->vehicle);

					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq21);
					addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
					myseqs.push_back(clusterU->seq12);
					myseqs.push_back(ysuccessor->seqi_n);
					resultMoves[2][3] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
			}

			if (!turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
				myseqs.push_back(clusterU->seq21);
				myseqs.push_back(clusterV->seqi_n);
				resultMoves[3][0] = seq->evaluationLB(myseqs,routeU->vehicle);
			}
			if (!clusterV->isDepot) // exchange X,U and V
			{
				if (!turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq1);
					addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
					myseqs.push_back(clusterU->seq21);
					myseqs.push_back(y->seqi_n);
					resultMoves[3][1] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
				if (!y->isDepot) // exchange X,U and (V,Y or Y,V)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq12);
					addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
					myseqs.push_back(clusterU->seq21);
					myseqs.push_back(ysuccessor->seqi_n);
					resultMoves[3][2] = seq->evaluationLB(myseqs,routeU->vehicle);

					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq21);
					addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
					myseqs.push_back(clusterU->seq21);
					myseqs.push_back(ysuccessor->seqi_n);
					resultMoves[3][3] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
			}
		}
	}
	else if (decalage == 2) // U and V are almost consecutive, we are in this situation : UXVY. Useful moves must be tested as special cases
	{
		// XUVY
		resultMoves[1][0] = seq->evaluationLB(Upred->seq0_i,clusterU->seq21,clusterV->seqi_n,routeU->vehicle);
		if (!clusterV->isDepot)
		{
			// VXUY
			resultMoves[1][1] = seq->evaluationLB(Upred->seq0_i,clusterV->seq1,clusterU->seq21,y->seqi_n,routeU->vehicle);
			// VUXY
			resultMoves[2][1] = seq->evaluationLB(Upred->seq0_i,clusterV->seq1,clusterU->seq12,y->seqi_n,routeU->vehicle);
			if (!y->isDepot)
			{
				// VYXU
				resultMoves[1][2] = seq->evaluationLB(Upred->seq0_i,clusterV->seq12,clusterU->seq21,ysuccessor->seqi_n,routeU->vehicle);
				// YVXU
				resultMoves[1][3] = seq->evaluationLB(Upred->seq0_i,clusterV->seq21,clusterU->seq21,ysuccessor->seqi_n,routeU->vehicle);
				// VYUX
				resultMoves[2][2] = seq->evaluationLB(Upred->seq0_i,clusterV->seq12,clusterU->seq12,ysuccessor->seqi_n,routeU->vehicle);
				// YVUX
				resultMoves[2][3] = seq->evaluationLB(Upred->seq0_i,clusterV->seq21,clusterU->seq12,ysuccessor->seqi_n,routeU->vehicle);
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// HERE STARTS THE SECOND PHASE AFTER PREPROCESSING
	// ONLY NON-FILTERED MOVES THAT HAVE A CHANCE TO BE IMPROVING
	/////////////////////////////////////////////////////////////

	for (int i=0 ; i<4 ; i++ )
	{
		for (int j=0 ; j<4 ; j++ )
		{
			shouldBeTested[i][j] = (resultMoves[i][j] < costZero) ;
		}
	}
	shouldBeTested[0][0] = true  ;
	resultMoves[0][0] = costZero ;


	// Same procedure as previously, but with "seq->evaluation" instead of "seq->evaluationLB"
	if (decalage >= 3)
	{
		if (!clusterV->isDepot && turned)
		{

			if (shouldBeTested[0][1])
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(clusterV->seq1);
				addSeqDataInPieces(clusterU,clusterV->place-1-clusterU->place);
				myseqs.push_back(y->seqi_n);
				resultMoves[0][1] = seq->evaluation(myseqs,routeU->vehicle);
			}

			if (!y->isDepot) // exchange U and (V,Y or Y,V)
			{
				if (shouldBeTested[0][2])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq12);
					addSeqDataInPieces(clusterU,clusterV->place-1-clusterU->place);
					myseqs.push_back(ysuccessor->seqi_n);
					resultMoves[0][2] = seq->evaluation(myseqs,routeU->vehicle);
				}

				if (shouldBeTested[0][3])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq21);
					addSeqDataInPieces(clusterU,clusterV->place-1-clusterU->place);
					myseqs.push_back(ysuccessor->seqi_n);
					resultMoves[0][3] = seq->evaluation(myseqs,routeU->vehicle);
				}
			}
		}

		if (!turned && shouldBeTested[1][0])
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			addSeqDataInPieces(x,clusterV->place-1-x->place);
			myseqs.push_back(clusterU->seq1);
			myseqs.push_back(clusterV->seqi_n);
			resultMoves[1][0] = seq->evaluation(myseqs,routeU->vehicle);
		}
		if (!clusterV->isDepot) // exchange U and V
		{
			if (shouldBeTested[1][1])
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(clusterV->seq1);
				addSeqDataInPieces(x,clusterV->place-1-x->place);
				myseqs.push_back(clusterU->seq1);
				myseqs.push_back(y->seqi_n);
				resultMoves[1][1] = seq->evaluation(myseqs,routeU->vehicle);
			}

			if (!y->isDepot && turned) // exchange U and (V,Y or Y,V)
			{
				if (shouldBeTested[1][2])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq12);
					addSeqDataInPieces(x,clusterV->place-1-x->place);
					myseqs.push_back(clusterU->seq1);
					myseqs.push_back(ysuccessor->seqi_n);
					resultMoves[1][2] = seq->evaluation(myseqs,routeU->vehicle);
				}

				if (shouldBeTested[1][3])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq21);
					addSeqDataInPieces(x,clusterV->place-1-x->place);
					myseqs.push_back(clusterU->seq1);
					myseqs.push_back(ysuccessor->seqi_n);
					resultMoves[1][3] = seq->evaluation(myseqs,routeU->vehicle);
				}
			}
		}

		if (!x->isDepot)
		{
			if (shouldBeTested[2][0] && !turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
				myseqs.push_back(clusterU->seq12);
				myseqs.push_back(clusterV->seqi_n);
				resultMoves[2][0] = seq->evaluation(myseqs,routeU->vehicle);
			}

			if (!clusterV->isDepot) // exchange U,X and V
			{
				if (shouldBeTested[2][1] && !turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq1);
					addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
					myseqs.push_back(clusterU->seq12);
					myseqs.push_back(y->seqi_n);
					resultMoves[2][1] = seq->evaluation(myseqs,routeU->vehicle);
				}
				if (!y->isDepot) // exchange U,X and (V,Y or Y,V)
				{
					if (shouldBeTested[2][2])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(clusterV->seq12);
						addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
						myseqs.push_back(clusterU->seq12);
						myseqs.push_back(ysuccessor->seqi_n);
						resultMoves[2][2] = seq->evaluation(myseqs,routeU->vehicle);
					}

					if (shouldBeTested[2][3])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(clusterV->seq21);
						addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
						myseqs.push_back(clusterU->seq12);
						myseqs.push_back(ysuccessor->seqi_n);
						resultMoves[2][3] = seq->evaluation(myseqs,routeU->vehicle);
					}
				}
			}

			if (shouldBeTested[3][0] && !turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
				myseqs.push_back(clusterU->seq21);
				myseqs.push_back(clusterV->seqi_n);
				resultMoves[3][0] = seq->evaluation(myseqs,routeU->vehicle);
			}
			if (!clusterV->isDepot) // exchange X,U and V
			{
				if (shouldBeTested[3][1] && !turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(clusterV->seq1);
					addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
					myseqs.push_back(clusterU->seq21);
					myseqs.push_back(y->seqi_n);
					resultMoves[3][1] = seq->evaluation(myseqs,routeU->vehicle);
				}
				if (!y->isDepot) // exchange X,U and (V,Y or Y,V)
				{
					if (shouldBeTested[3][2])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(clusterV->seq12);
						addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
						myseqs.push_back(clusterU->seq21);
						myseqs.push_back(ysuccessor->seqi_n);
						resultMoves[3][2] = seq->evaluation(myseqs,routeU->vehicle);
					}

					if (shouldBeTested[3][3])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(clusterV->seq21);
						addSeqDataInPieces(xsuccessor,clusterV->place-1-xsuccessor->place);
						myseqs.push_back(clusterU->seq21);
						myseqs.push_back(ysuccessor->seqi_n);
						resultMoves[3][3] = seq->evaluation(myseqs,routeU->vehicle);
					}
				}
			}
		}
	}
	else if (decalage == 2)
	{
		// XUVY
		if (shouldBeTested[1][0]) resultMoves[1][0] = seq->evaluation(Upred->seq0_i,clusterU->seq21,clusterV->seqi_n,routeU->vehicle);
		if (!clusterV->isDepot)
		{
			// VXUY
			if (shouldBeTested[1][1]) resultMoves[1][1] = seq->evaluation(Upred->seq0_i,clusterV->seq1,clusterU->seq21,y->seqi_n,routeU->vehicle);
			// VUXY
			if (shouldBeTested[2][1]) resultMoves[2][1] = seq->evaluation(Upred->seq0_i,clusterV->seq1,clusterU->seq12,y->seqi_n,routeU->vehicle);
			if (!y->isDepot)
			{
				// VYXU
				if (shouldBeTested[1][2]) resultMoves[1][2] = seq->evaluation(Upred->seq0_i,clusterV->seq12,clusterU->seq21,ysuccessor->seqi_n,routeU->vehicle);
				// YVXU
				if (shouldBeTested[1][3]) resultMoves[1][3] = seq->evaluation(Upred->seq0_i,clusterV->seq21,clusterU->seq21,ysuccessor->seqi_n,routeU->vehicle);
				// VYUX
				if (shouldBeTested[2][2]) resultMoves[2][2] = seq->evaluation(Upred->seq0_i,clusterV->seq12,clusterU->seq12,ysuccessor->seqi_n,routeU->vehicle);
				// YVUX
				if (shouldBeTested[2][3]) resultMoves[2][3] = seq->evaluation(Upred->seq0_i,clusterV->seq21,clusterU->seq12,ysuccessor->seqi_n,routeU->vehicle);
			}
		}
	}

	// Identify the best move involving U,V and apply
	ibest = 0 ;
	jbest = 0 ;
	moveMin = 1.e30 ;
	for (int i=0 ; i<4 ; i++)
	{
		for (int j=0 ; j<4 ; j++ )
		{
			if (shouldBeTested[i][j] && resultMoves[i][j] < moveMin - EPSILON_LS )
			{
				moveMin = resultMoves[i][j] ;
				ibest = i ;
				jbest = j ;
			}
		}
	}

	// No improving move has been found
	if ( ibest == 0 && jbest == 0 ) 
	{ 
		clusterU = tempU ;
		clusterV = tempV ;
		x = clusterU->successor ;
		y = clusterV->successor ;
		return 0 ;
	}


	///////////////////////////////////////////
	// HERE STARTS THE MOVE APPLICATION PROCESS
	///////////////////////////////////////////

	reinitMoves(routeU);  // Say that all moves involving this route must be tested again

	// Update the solution
	if (decalage >= 3)
	{
		if ( ibest == 1 || ibest == 2) insertClusterInRoute(clusterU,placeV);
		if ( ibest == 2 ) insertClusterInRoute(x,clusterU);
		if ( ibest == 3 ) { insertClusterInRoute(x,placeV); insertClusterInRoute(clusterU,x); }

		if ( jbest == 1 || jbest == 2) insertClusterInRoute(clusterV,placeU);
		if ( jbest == 2) insertClusterInRoute(y,clusterV);
		if ( jbest == 3 ) { insertClusterInRoute(y,placeU); insertClusterInRoute(clusterV,y); }
	}

	// Special cases of decalage == 2
	else if (decalage == 2)
	{
		if (ibest == 1 && jbest == 0) { insertClusterInRoute(clusterU,x); }
		else if (ibest == 1 && jbest == 1) { insertClusterInRoute(x,clusterV); insertClusterInRoute(clusterU,x);  }
		else if (ibest == 1 && jbest == 2) { insertClusterInRoute(x,y); insertClusterInRoute(clusterU,x);  }
		else if (ibest == 1 && jbest == 3) { insertClusterInRoute(clusterV,y); insertClusterInRoute(x,clusterV); insertClusterInRoute(clusterU,x);  }
		else if (ibest == 2 && jbest == 1) { insertClusterInRoute(clusterV,clusterU->pred) ;  }
		else if (ibest == 2 && jbest == 2) { insertClusterInRoute(clusterV,clusterU->pred) ;  insertClusterInRoute(y,clusterV) ; }
		else if (ibest == 2 && jbest == 3) { insertClusterInRoute(y,clusterU->pred) ;  insertClusterInRoute(clusterV,y) ; }
		else throw string ("ERROR move intra-route") ;
	} 

	routeU->updateRouteData(false); // Update the pre-processed data on the subsequences of the route
	setRouteVide(); // Keep a pointer on the first empty route
	rechercheTerminee = false ; // Not finished the search
	nbIntraSwap ++ ;
	clusterU = tempU ;
	clusterV = tempV ;
	x = clusterU->successor ;
	y = clusterV->successor ;
	return 1 ; // Return Success
}


int LocalSearch::intraRoute2Opt ()
{
	// Evaluation procedure for 2-Opt
	ClusterInRoute * nodeNum = clusterU->successor ;
	ClusterInRoute * nodeUpred = clusterU->pred ;
	ClusterInRoute * temp ;
	SeqData * seq = clusterU->seq0_i ;

	double cost ;
	double costZero = routeU->currentRouteCost ;

	myseqs.clear();
	myseqs.push_back(clusterU->pred->seq0_i);
	addReverseSeqDataInPieces(clusterU,clusterV->place-clusterU->place);
	myseqs.push_back(clusterV->successor->seqi_n);

	// Compute the lower bound on move value and exits if no possible improvement
	cost = seq->evaluationLB(myseqs,routeU->vehicle) ;
	if (cost - costZero  > -EPSILON_LS)  
		return 0 ;

	// Compute the real move value
	cost = seq->evaluation(myseqs,routeU->vehicle) ;
	if (cost - costZero  > -EPSILON_LS)  
		return 0 ;

	// Apply the move and update the solution
	reinitMoves(routeU);
	clusterU->pred = nodeNum ;
	clusterU->successor = y ;

	while ( nodeNum != clusterV )
	{
		temp = nodeNum->successor ;
		nodeNum->successor = nodeNum->pred ;
		nodeNum->pred = temp ;
		nodeNum = temp ;
	}

	clusterV->successor = clusterV->pred ;
	clusterV->pred = nodeUpred ;
	nodeUpred->successor = clusterV ;
	y->pred = clusterU ;
	routeU->updateRouteData(false);  // Update the pre-processed data on the subsequences of the route
	setRouteVide(); // Keep a pointer on the first empty route
	rechercheTerminee = false ; // Not finished the search
	nbIntra2Opt ++ ;
	return 1 ; // Return Success
}

void LocalSearch::shuffleDepotAndCustomerIndices ()
{
	// Shuffling the shuffledDepotAndCustomerIndices vector
	int j,temp ;
  for (int i = 0 ; i < (int)shuffledDepotAndCustomerIndices.size() - 1 ; i++)
  {
    j = i + rand() % ((int)shuffledDepotAndCustomerIndices.size() - i) ;
    temp = shuffledDepotAndCustomerIndices[i] ;
    shuffledDepotAndCustomerIndices[i] = shuffledDepotAndCustomerIndices[j] ;
    shuffledDepotAndCustomerIndices[j] = temp ;
  }
}

void LocalSearch::updateMoves ()
{
	int client, client2, size ;
  for (int i=0 ; i<( int)shuffledDepotAndCustomerIndices.size() ; i++)
  {
    client = shuffledDepotAndCustomerIndices[i] ;
    clients[client].moves.clear();
    size = (int)params->cli[client].inverseCloseCustomers.size();
    for (int a1 = 0 ; a1 < size ; a1++ )
    {
      client2 = params->cli[client].inverseCloseCustomers[a1] ;
      if (client2 >= params->nbDepots && clients[client2].isPresent) clients[client].moves.push_back(client2);
    }
  }
}

void LocalSearch::setRouteVide()
{
	emptyRoute = NULL ;
	int route = 0 ;
	while (emptyRoute == NULL && route < params->nbVehiculesPerDep )
	{
		if (depots[route].successor->isDepot) emptyRoute = depots[route].route ;
		route ++ ;
	}
}

void LocalSearch::nodeTestedForEachRoute (int cli)
{
	for (int route = 0 ; route < (int)params->nombreVehicules ; route ++)
		routes[route].nodeAndRouteTested[cli]=true ;	
}

void LocalSearch::insertClusterInRoute(ClusterInRoute * U, ClusterInRoute * V)
{
	if (U->pred != V && U != V)
	{
		U->pred->successor = U->successor ;
		U->successor->pred = U->pred ;
		V->successor->pred = U ;
		U->pred = V ;
		U->successor = V->successor ;
		V->successor = U ;
		U->route = V->route ;
	}
}

void LocalSearch::reinitMoves(Route * r)
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
		r->nodeAndRouteTested[i] = false ;

	for (ClusterInRoute * tempClusterInRoute = r->depot->successor ; !tempClusterInRoute->isDepot ; tempClusterInRoute = tempClusterInRoute->successor)
		for (int route = 0 ; route < params->nbVehiculesPerDep ; route ++)
			depots[route].route->nodeAndRouteTested[tempClusterInRoute->cluster_idx] = false ;
}

void LocalSearch::reinitAllMoves()
{
  for (int route = 0 ; route < params->nbVehiculesPerDep ; route ++)
    for (int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
      depots[route].route->nodeAndRouteTested[i] = false ;
}

void LocalSearch::addSeqDataInPieces (ClusterInRoute * node, int length)
{
	ClusterInRoute * curCluster = node ;
	SeqData * curSeq ;
	int curLength = length + 1 ;
	int temp ;
	while (curLength > 0)
	{
		temp = min(curLength,params->sizeSD) ;
		curSeq = curCluster->seqi_j[temp-1];
		myseqs.push_back(curSeq);
		curLength -= temp ;
		if (curLength != 0) curCluster = clients[curSeq->lastCluster].successor ;
	}
}

void LocalSearch::addReverseSeqDataInPieces (ClusterInRoute * node, int length)
{
	ClusterInRoute * curCluster = node ;
	SeqData * curSeq ;
	list <SeqData*> myseqTemp ;
	int curLength = length + 1 ;
	int temp ;
	while (curLength > 0)
	{
		temp = min(curLength,params->sizeSD) ;
		curSeq = curCluster->seqj_i[temp-1];
		myseqTemp.push_back(curSeq);
		curLength -= temp ;
		if (curLength != 0) curCluster = clients[curSeq->firstCluster].successor ;
	}

	for ( list<SeqData*>::reverse_iterator it=myseqTemp.rbegin() ; it != myseqTemp.rend() ; it++)
		myseqs.push_back(*it);
}

LocalSearch::LocalSearch(void)
{
	allAttributesSet = false ;
	seqdeb = NULL ;
}

LocalSearch::LocalSearch(Params * params,Individu * individu) : params (params),individu(individu)
{
	int nbVeh ;
	allAttributesSet = true ;
	nbInterSwap = 0 ;
	nbIntraSwap = 0 ;
	nbInter2Opt = 0 ;
	nbIntra2Opt = 0 ;
	seqdeb = NULL ;
	nbTotalRISinceBeginning = 0 ;
	vector < ClusterInRoute * > tempClusterInRoute ;
	vector <int> temp2 ;

	nbVeh = params->nombreVehicules ;
	clients = new ClusterInRoute [params->nbClients + params->nbDepots + 1] ;
	depots = new ClusterInRoute [nbVeh] ;
	depotsFin = new ClusterInRoute [nbVeh] ;
	routes = new Route [nbVeh] ; 

	for (int i = 0 ; i <  params->nbClients + params->nbDepots ; i ++ )
		clients[i] = ClusterInRoute(false,i,false,NULL,NULL,NULL,params);

	for (int i = 0 ; i < nbVeh ; i++ )
	{
		depots[i] = ClusterInRoute(true,params->ordreVehicules[i].depotNumber,false,NULL,NULL,NULL,params);
		depotsFin[i] = ClusterInRoute(true,params->ordreVehicules[i].depotNumber,false,NULL,NULL,NULL,params);
		routes[i] = Route(i,0,&params->ordreVehicules[i],params,individu);
		depots[i].route = &routes[i] ;
		depotsFin[i].route = &routes[i] ;
		routes[i].depot = &depots[i] ;
	}
	

	for (int i=0 ; i < 4 ; i++)
	{
		resultMoves.push_back(vector<double>(4));
		shouldBeTested.push_back(vector<bool>(4));
	}

  emptyRoute = NULL;

	for (int i=params->nbDepots ; i < params->nbDepots + params->nbClients ; i++)
		shuffledDepotAndCustomerIndices.push_back(i);

	// Initialization of the SeqData structures, for all nodes
	// These structures were designed to stay contiguous in memory (however it had only a negligible impact on performance)
	int nbSeqsSet = 0 ;
	int taillemyseqDatas = (params->sizeSD + params->sizeSD + 4)*(params->nbClients + params->nbDepots + 2*params->nbVehiculesPerDep) ;
	SeqData * myseqDatas = new SeqData [taillemyseqDatas] ;
	seqdeb = myseqDatas ;
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i++) 
	{
		for (int ii=nbSeqsSet ; ii < nbSeqsSet+4+params->sizeSD+params->sizeSD ; ii++)
			myseqDatas[ii].initialisation(clients[i].cluster_idx,params,individu,false);  

		clients[i].seq0_i = &myseqDatas[nbSeqsSet] ;
		clients[i].seqi_n = &myseqDatas[nbSeqsSet+1] ;
		clients[i].seqi_0 = &myseqDatas[nbSeqsSet+2] ;
		clients[i].seqn_i = &myseqDatas[nbSeqsSet+3] ;
		for (int ii=0 ; ii<params->sizeSD ; ii++)
			clients[i].seqi_j.push_back(&myseqDatas[nbSeqsSet+4+ii]);
		for (int ii=0 ; ii<params->sizeSD ; ii++)
			clients[i].seqj_i.push_back(&myseqDatas[nbSeqsSet+4+params->sizeSD+ii]);
		clients[i].setRemaining();
		nbSeqsSet += 4+params->sizeSD + params->sizeSD ;
	}

	for (int i=0 ; i < params->nbVehiculesPerDep  ; i++) 
	{
		for (int ii=nbSeqsSet ; ii < nbSeqsSet+4+params->sizeSD+params->sizeSD ; ii++)
			myseqDatas[ii].initialisation(depots[i].cluster_idx,params,individu,false);

		depots[i].seq0_i = &myseqDatas[nbSeqsSet] ;
		depots[i].seqi_n = &myseqDatas[nbSeqsSet+1] ;
		depots[i].seqi_0 = &myseqDatas[nbSeqsSet+2] ;
		depots[i].seqn_i = &myseqDatas[nbSeqsSet+3] ;
		for (int ii=0 ; ii<params->sizeSD ; ii++)
			depots[i].seqi_j.push_back(&myseqDatas[nbSeqsSet+4+ii]);
		for (int ii=0 ; ii<params->sizeSD ; ii++)
			depots[i].seqj_i.push_back(&myseqDatas[nbSeqsSet+4+params->sizeSD+ii]);
		depots[i].setRemaining();
		nbSeqsSet += 4+params->sizeSD + params->sizeSD ;
	}

	for (int i=0 ; i < params->nbVehiculesPerDep ; i++) 
	{
		for (int ii=nbSeqsSet ; ii < nbSeqsSet+4+params->sizeSD+params->sizeSD ; ii++)
			myseqDatas[ii].initialisation(depotsFin[i].cluster_idx,params,individu,false);

		depotsFin[i].seq0_i = &myseqDatas[nbSeqsSet] ;
		depotsFin[i].seqi_n = &myseqDatas[nbSeqsSet+1] ;
		depotsFin[i].seqi_0 = &myseqDatas[nbSeqsSet+2] ;
		depotsFin[i].seqn_i = &myseqDatas[nbSeqsSet+3] ;
		for (int ii=0 ; ii<params->sizeSD ; ii++)
			depotsFin[i].seqi_j.push_back(&myseqDatas[nbSeqsSet+4+ii]);
		for (int ii=0 ; ii<params->sizeSD ; ii++)
			depotsFin[i].seqj_i.push_back(&myseqDatas[nbSeqsSet+4+params->sizeSD+ii]);
		depotsFin[i].setRemaining();
		nbSeqsSet += 4+params->sizeSD + params->sizeSD ;
	}
	
}

LocalSearch::~LocalSearch(void)
{	
	SeqData * seqdeb2 = (SeqData*) seqdeb ;
	delete [] seqdeb2 ;

	if (allAttributesSet)
	{
		delete [] clients ;
		delete [] depots ;
		delete [] depotsFin ;
		delete [] routes ;
	}
}

