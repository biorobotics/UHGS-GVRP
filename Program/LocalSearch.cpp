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
	melangeParcours();
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
	Noeud * tempNoeud ;
	int size2 ;

	// We search and apply moves until a local minimum is attained
	while (!rechercheTerminee)
	{
		rechercheTerminee = true ;
		moveEffectue = 0 ;

		// For every node U in random order
		for ( int posU = 0 ; posU < (int)ordreParcours.size() ; posU ++ )
		{
			posU -= moveEffectue ; // We return on the last node if there has been a move
			nbMoves += moveEffectue ;
			moveEffectue = 0 ;
			routeVideTestee = false ;

			noeudU = &clients[ordreParcours[posU]];
			routeU = noeudU->route ;
			vehicleU = routeU->vehicle ;
			x = noeudU->suiv ;

			// In the CARP, some service removal do not reduce the cost of the route
			// In this case, very few moves can improve the solution (only 2-opt variants), and thus SWAP and RELOCATE variants do not need to be tested
			costBeforeRemoval = noeudU->seq0_i->evaluation(routeU->depot->pred->seq0_i,routeU->vehicle);
			costAfterRemoval = noeudU->seq0_i->evaluation(noeudU->pred->seq0_i,noeudU->suiv->seqi_n,routeU->vehicle);
			gainWhenRemoving = (costAfterRemoval < costBeforeRemoval - 0.1) ;

			// For every node V in random order
			size2 = (int)noeudU->moves.size();
			for ( int posV = 0 ; posV < size2 && moveEffectue == 0 ; posV ++ )
			{
				noeudV = &clients[noeudU->moves[posV]] ;
				routeV = noeudV->route ;
				vehicleV = routeV->vehicle ;

				// If we have not yet tested the moves involving the node U and the route of node V
				// (This flag is reset to false as soon as there is a modification in the route)
				if (!noeudV->route->nodeAndRouteTested[noeudU->customer_or_depot_idx])
				{
					y = noeudV->suiv ;
					if (routeV->customer_or_depot_idx != routeU->customer_or_depot_idx)
					{
						if (moveEffectue != 1)
						{
							tempNoeud = noeudV ;
							noeudV = noeudV->suiv ;
							y = noeudV->suiv ;
							// Testing Relocate, Swap, CROSS and I-CROSS (limited to 2 customers) of nodeU and nodeV 
							// Case where they are in different routes
							if (gainWhenRemoving) moveEffectue = interRouteGeneralInsert(); 
							noeudV = tempNoeud ;
							y = noeudV->suiv ;
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
						tempNoeud = noeudV ;
						noeudV = noeudV->suiv ;
						y = noeudV->suiv ;

						// Testing Relocate, Swap, CROSS and I-CROSS (limited to 2 customers) of nodeU and nodeV 
						// Case where they are in the same route
						if (moveEffectue != 1 && gainWhenRemoving) 
							moveEffectue = intraRouteGeneralInsertDroite();
						noeudV = tempNoeud ;
						y = noeudV->suiv ;
					}
				}
			}

			noeudV = noeudU->suiv ;
			routeV = noeudV->route ;
			vehicleV = routeV->vehicle ;
			y = noeudV->suiv ;
			if (!noeudV->route->nodeAndRouteTested[noeudU->customer_or_depot_idx])
			{
				while (moveEffectue != 1 && !noeudV->estUnDepot)
				{ 
					// Testing 2-Opt between U and V (if the restriction of the granular search allows) 
					if (params->isCorrelated[noeudU->pred->customer_or_depot_idx][noeudV->customer_or_depot_idx] || params->isCorrelated[noeudU->customer_or_depot_idx][noeudV->suiv->customer_or_depot_idx]) 
						moveEffectue = intraRoute2Opt ();
					noeudV = noeudV->suiv ;
					y = noeudV->suiv ;
				}
			}

			// Special cases : testing the insertions behind the depot, and the empty routes
			for (int route = 0 ; route < params->nbVehiculesPerDep && moveEffectue == 0 ; route ++)
			{
				noeudV = &depots[route] ;
				routeV = noeudV->route ;
				y = noeudV->suiv ;
				if ( (!noeudV->route->nodeAndRouteTested[noeudU->customer_or_depot_idx]) && (!y->estUnDepot || !routeVideTestee))
				{
					if (y->estUnDepot) routeVideTestee = true ;
					if (routeV != routeU)
					{
						tempNoeud = noeudV ;
						noeudV = depots[route].suiv ;
						y = noeudV->suiv ;

						// Insertion after the depot, in a different route
						if (gainWhenRemoving && (params->isCorrelated[noeudU->customer_or_depot_idx][noeudV->customer_or_depot_idx] || params->isCorrelated[noeudU->customer_or_depot_idx][y->customer_or_depot_idx]) && moveEffectue != 1 ) 
							moveEffectue = interRouteGeneralInsert();

						noeudV = noeudV->route->depot ;
						y = noeudV->suiv ;

						// 2-Opt* after the depot
						if (params->isCorrelated[noeudU->pred->customer_or_depot_idx][noeudV->customer_or_depot_idx] && moveEffectue != 1) 
							moveEffectue = interRoute2Opt ();

						// 2-Opt* after the depot
						if ((params->isCorrelated[x->customer_or_depot_idx][y->customer_or_depot_idx] || params->isCorrelated[y->customer_or_depot_idx][x->customer_or_depot_idx]) && moveEffectue != 1) 
							moveEffectue = interRoute2OptInv ();

						noeudV = tempNoeud ;
						y = noeudV->suiv ;
					}
					else
					{
						tempNoeud = noeudV ;
						noeudV = depots[route].suiv ;
						y = noeudV->suiv ;

						// Insertion after the depot, in the same route
						if ((params->isCorrelated[noeudU->customer_or_depot_idx][noeudV->customer_or_depot_idx] || params->isCorrelated[noeudU->customer_or_depot_idx][y->customer_or_depot_idx]) && moveEffectue != 1) 
							moveEffectue = intraRouteGeneralInsertDroite();

						noeudV = tempNoeud ;
						y = noeudV->suiv ;
					}
				}
			}

			// Say that we have tested the node U with all routes
			if (moveEffectue == 0)
				nodeTestedForEachRoute(noeudU->customer_or_depot_idx);
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
	SeqData * seq = noeudU->seq0_i ;
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
	if (!noeudU->estUnDepot)
	{
		// We evaluate the result of the move as a combination of existing subsequences
		resultMoves[1][0] = seq->evaluationLB(noeudU->pred->seq0_i,x->seqi_n,routeU->vehicle);

		// Can we receive something (if V is not a depot)
		// In the following, all the conditionals are set up to avoid moving depots
		if (!noeudV->estUnDepot)
			resultMoves[1][1] = seq->evaluationLB(noeudU->pred->seq0_i,noeudV->seq1,x->seqi_n,routeU->vehicle);

		// And so on...
		if (!x->estUnDepot)
		{
			resultMoves[2][0] = seq->evaluationLB(noeudU->pred->seq0_i,x->suiv->seqi_n,routeU->vehicle);
			resultMoves[3][0] = resultMoves[2][0] ;
			if (!noeudV->estUnDepot)
			{
				resultMoves[2][1] = seq->evaluationLB(noeudU->pred->seq0_i,noeudV->seq1,x->suiv->seqi_n,routeU->vehicle);
				resultMoves[3][1] = resultMoves[2][1] ;
				if (!y->estUnDepot)
				{
					resultMoves[2][2] = seq->evaluationLB(noeudU->pred->seq0_i,noeudV->seq12,x->suiv->seqi_n,routeU->vehicle);
					resultMoves[2][3] = seq->evaluationLB(noeudU->pred->seq0_i,noeudV->seq21,x->suiv->seqi_n,routeU->vehicle);
					resultMoves[3][2] = resultMoves[2][2];
					resultMoves[3][3] = resultMoves[2][3];
				}
			}
		}
	}

	if (!noeudU->estUnDepot)
	{
		resultMoves[1][0] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq1,noeudV->seqi_n,routeV->vehicle);
		if (!x->estUnDepot)
		{
			resultMoves[2][0] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq12,noeudV->seqi_n,routeV->vehicle);
			resultMoves[3][0] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq21,noeudV->seqi_n,routeV->vehicle);
		}
	}

	if (!noeudV->estUnDepot)
	{
		if (!noeudU->estUnDepot)
		{
			resultMoves[1][1] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq1,y->seqi_n,routeV->vehicle);
			if (!x->estUnDepot)
			{
				resultMoves[2][1] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq12,y->seqi_n,routeV->vehicle);
				resultMoves[3][1] += seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq21,y->seqi_n,routeV->vehicle);
			}
		}
		if (!y->estUnDepot && !noeudU->estUnDepot && !x->estUnDepot)
		{
			temp = seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq12,y->suiv->seqi_n,routeV->vehicle);
			resultMoves[2][2] += temp ;
			resultMoves[2][3] += temp ;
			temp = seq->evaluationLB(noeudV->pred->seq0_i,noeudU->seq21,y->suiv->seqi_n,routeV->vehicle);
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

	if (!noeudU->estUnDepot)
	{
		if (shouldBeTested[1][0]) resultMoves[1][0] = seq->evaluation(noeudU->pred->seq0_i,x->seqi_n,routeU->vehicle);

		if (!noeudV->estUnDepot && shouldBeTested[1][1])
			resultMoves[1][1] = seq->evaluation(noeudU->pred->seq0_i,noeudV->seq1,x->seqi_n,routeU->vehicle);

		if (!x->estUnDepot)
		{
			if (shouldBeTested[2][0] || shouldBeTested[3][0]) 
			{
				resultMoves[2][0] = seq->evaluation(noeudU->pred->seq0_i,x->suiv->seqi_n,routeU->vehicle);
				resultMoves[3][0] = resultMoves[2][0] ;
			}
			if (!noeudV->estUnDepot)
			{
				if (shouldBeTested[2][1] || shouldBeTested[3][1])
				{
					resultMoves[2][1] = seq->evaluation(noeudU->pred->seq0_i,noeudV->seq1,x->suiv->seqi_n,routeU->vehicle);
					resultMoves[3][1] = resultMoves[2][1] ;
				}
				if (!y->estUnDepot)
				{
					if (shouldBeTested[2][2] || shouldBeTested[3][2])
					{
						resultMoves[2][2] = seq->evaluation(noeudU->pred->seq0_i,noeudV->seq12,x->suiv->seqi_n,routeU->vehicle);
						resultMoves[3][2] = resultMoves[2][2];
					}
					if (shouldBeTested[2][3] || shouldBeTested[3][3])
					{
						resultMoves[2][3] = seq->evaluation(noeudU->pred->seq0_i,noeudV->seq21,x->suiv->seqi_n,routeU->vehicle);
						resultMoves[3][3] = resultMoves[2][3];
					}
				}
			}
		}
	}

	if (!noeudU->estUnDepot)
	{
		if (shouldBeTested[1][0]) resultMoves[1][0] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq1,noeudV->seqi_n,routeV->vehicle);
		if (!x->estUnDepot)
		{
			if (shouldBeTested[2][0]) resultMoves[2][0] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq12,noeudV->seqi_n,routeV->vehicle);
			if (shouldBeTested[3][0]) resultMoves[3][0] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq21,noeudV->seqi_n,routeV->vehicle);
		}
	}

	if (!noeudV->estUnDepot)
	{
		if (!noeudU->estUnDepot)
		{
			if (shouldBeTested[1][1]) resultMoves[1][1] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq1,y->seqi_n,routeV->vehicle);
			if (!x->estUnDepot)
			{
				if (shouldBeTested[2][1]) resultMoves[2][1] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq12,y->seqi_n,routeV->vehicle);
				if (shouldBeTested[3][1]) resultMoves[3][1] += seq->evaluation(noeudV->pred->seq0_i,noeudU->seq21,y->seqi_n,routeV->vehicle);
			}
		}
		if (!y->estUnDepot && !noeudU->estUnDepot && !x->estUnDepot)
		{
			if (shouldBeTested[2][2] || shouldBeTested[2][3])
			{
				temp = seq->evaluation(noeudV->pred->seq0_i,noeudU->seq12,y->suiv->seqi_n,routeV->vehicle);
				resultMoves[2][2] += temp ;
				resultMoves[2][3] += temp ;
			}
			if (shouldBeTested[3][2] || shouldBeTested[3][3])
			{
				temp = seq->evaluation(noeudV->pred->seq0_i,noeudU->seq21,y->suiv->seqi_n,routeV->vehicle);
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

	Noeud * placeV = noeudV->pred ;
	Noeud * placeU = noeudU->pred ;
	reinitMoves(placeU->route); // Say that all moves involving this route must be tested again
	reinitMoves(placeV->route); // Say that all moves involving this route must be tested again

	// Moving the nodes in the solution structure
	if ( ibest == 1 || ibest == 2) insertNoeud(noeudU,placeV);
	if ( ibest == 2) insertNoeud(x,noeudU);
	if ( ibest == 3 ) { insertNoeud(x,placeV); insertNoeud(noeudU,x); }

	if ( jbest == 1 || jbest == 2) insertNoeud(noeudV,placeU);
	if ( jbest == 2) insertNoeud(y,noeudV);
	if ( jbest == 3 ) { insertNoeud(y,placeU); insertNoeud(noeudV,y); }

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
	SeqData * seq = noeudU->seq0_i ;

	if  (routeU->depot->customer_or_depot_idx != routeV->depot->customer_or_depot_idx || (noeudU->pred->estUnDepot && noeudV->estUnDepot)) 
		return 0 ; // Cannot do a 2-Opt* if the depot is placed in the wrong way

	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ; // Cost before the move

	// Lower bound on the move value
	cost = seq->evaluationLB(noeudU->pred->seq0_i,y->seqi_n,routeU->vehicle) + seq->evaluationLB(noeudV->seq0_i,noeudU->seqi_n,routeV->vehicle) - costZero ;
	if ( cost  > -EPSILON_LS ) 
		return 0 ; // Exit if no chance of improvement

	// Exact move evaluation
	cost = seq->evaluation(noeudU->pred->seq0_i,y->seqi_n,routeU->vehicle) + seq->evaluation(noeudV->seq0_i,noeudU->seqi_n,routeV->vehicle) - costZero ;
	if ( cost  > -EPSILON_LS ) 
		return 0 ; 

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	reinitMoves(routeU);  // Say that all moves involving this route must be tested again
	reinitMoves(routeV);  // Say that all moves involving this route must be tested again
	Noeud * tempU = noeudU ;
	noeudU = noeudU->pred ;
	x = noeudU->suiv ;

	// Updating the solution
	Noeud * count ;
	Noeud * depotU = routeU->depot ;
	Noeud * depotV = routeV->depot ;
	Noeud * depotUFin = depotU->pred ;
	Noeud * depotVFin = depotV->pred ;
	Noeud * depotUpred = depotUFin->pred ;

	count = y ;
	while ( !count->estUnDepot )
	{
		count->route = routeU ;
		count = count->suiv ;
	}

	count = x ;
	while ( !count->estUnDepot )
	{
		count->route = routeV ;
		count = count->suiv ;
	}

	noeudU->suiv = y ;
	y->pred = noeudU ;
	noeudV->suiv = x ;
	x->pred = noeudV ;

	if (x->estUnDepot)
	{
		depotUFin->pred = depotVFin->pred ;
		depotUFin->pred->suiv = depotUFin ;
		noeudV->suiv = depotVFin ;
		depotVFin->pred = noeudV ;
	}
	else
	{
		depotUFin->pred = depotVFin->pred ;
		depotUFin->pred->suiv = depotUFin ;
		depotVFin->pred = depotUpred ;
		depotVFin->pred->suiv = depotVFin ;
	}

	// Update the pre-processed data on the subsequences of the route
	routeU->updateRouteData(false);
	routeV->updateRouteData(false);
	setRouteVide(); // Keep a pointer on the first empty route

	rechercheTerminee = false ; // Not finished the search
	nbInter2Opt ++ ;
	noeudU = tempU ;
	x = noeudU->suiv ;
	return 1 ; // Return Success
}


int LocalSearch::interRoute2OptInv()
{
	// 2-Opt* with route inversions
	SeqData * seq = noeudU->seq0_i ;
	double cost ;
	double costTemp ;
	double costTempReverse = 1.e30 ;
	bool reverseRouteU, reverseRouteV ; 
	if (noeudU->route == noeudV->route) { return 0 ; }
	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ;

	// Compute a lower bound on the value of the move
	cost = - costZero ;
	costTemp = seq->evaluationLB(y->seqn_i,x->seqi_n,routeU->vehicle) ;
	costTempReverse = seq->evaluationLB(x->seqn_i,y->seqi_n,routeU->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteU = false ; }
	else  { cost += costTempReverse ; reverseRouteU = true ; }

	costTemp = seq->evaluationLB(noeudV->seq0_i,noeudU->seqi_0,routeV->vehicle) ;
	costTempReverse = seq->evaluationLB(noeudU->seq0_i,noeudV->seqi_0,routeV->vehicle) ;
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

	costTemp = seq->evaluation(noeudV->seq0_i,noeudU->seqi_0,routeV->vehicle) ;
	costTempReverse = seq->evaluation(noeudU->seq0_i,noeudV->seqi_0,routeV->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteV = false ; }
	else  { cost += costTempReverse ; reverseRouteV = true ; }

	if ( cost  > -EPSILON_LS ) 
		return 0 ;

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	reinitMoves(noeudU->route);  // Say that all moves involving this route must be tested again
	reinitMoves(noeudV->route);  // Say that all moves involving this route must be tested again
	Noeud * depotU = routeU->depot ;
	Noeud * depotV = routeV->depot ;
	Noeud * depotUFin = routeU->depot->pred ;
	Noeud * depotVFin = routeV->depot->pred ;
	Noeud * depotUSuiv = depotU->suiv ;

	// Update the solution
	Noeud * temp ;
	Noeud * yy = y ;
	Noeud * uu = noeudU ;

	while ( !yy->estUnDepot )
	{
		temp = yy->suiv ;
		yy->suiv = yy->pred ;
		yy->pred = temp ;
		yy->route = routeU ;
		yy = temp ;
	}

	while ( !uu->estUnDepot )
	{
		temp = uu->pred ;
		uu->pred = uu->suiv ;
		uu->suiv = temp ;
		uu->route = routeV ;
		uu = temp ;
	}

	noeudV->suiv = noeudU ;
	noeudU->pred = noeudV ;
	y->suiv = x ;
	x->pred = y ;

	if (y->estUnDepot)
	{
		depotVFin->suiv = depotV ;
		depotVFin->pred = depotUSuiv ;
		depotVFin->pred->suiv = depotVFin ;
		depotU->suiv = x ;
		x->pred = depotU ;
	}
	else if (noeudU->estUnDepot )
	{
		depotU->suiv = depotVFin->pred ;
		depotU->suiv->pred = depotU ;
		depotU->pred = depotUFin ;
		depotVFin->pred = noeudV ;
		noeudV->suiv = depotVFin ;
	}
	else
	{
		depotU->suiv = depotVFin->pred ;
		depotU->suiv->pred = depotU ;
		depotVFin->pred = depotUSuiv ;
		depotVFin->pred->suiv = depotVFin ;
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
	Noeud * tempU = noeudU ;
	Noeud * tempV = noeudV ;
	bool turned = false ;

	int decalage = noeudV->place - noeudU->place ; // Decalage is the difference of index between U and V
	if (decalage >= -1 && decalage <= 1) return 0 ; // This means that U and V are already consecutive, testing these moves is useless

	// If decalage < 0, this means that V is on the left of U in the route
	// We don't want to deal with this case, so we simply inverse the notations of U and V 
	// And we deal in the following with only the case where V is on the right
	if ( decalage < 0 )
	{
		noeudU = tempV ;
		noeudV = tempU ;
		x = noeudU->suiv ;
		y = noeudV->suiv ;
		decalage = -decalage ;
		turned = true ;
	}

	double moveMin ;
	int ibest, jbest ;
	SeqData * seq = noeudU->seq0_i ;
	Noeud * Upred = noeudU->pred ;
	Noeud * placeU = noeudU->pred ;
	Noeud * placeV = noeudV->pred ;
	Noeud * xsuiv = x->suiv ;
	Noeud * ysuiv = y->suiv ;
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
		if (!noeudV->estUnDepot && turned) // relocate of V before U (V cannot be a depot, but U can be)
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			myseqs.push_back(noeudV->seq1);
			// The intra-route moves are the only moves which can involve general subsequences (i,j) in the middle of the route
			// To reduce a bit the pre-processing, we keep only sequences of limited size (<= 15), even if this involves to concatenate more pieces
			// The concatenation of the good number of pieces is done by "addSeqDataInPieces"
			addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place);
			myseqs.push_back(y->seqi_n);
			resultMoves[0][1] = seq->evaluationLB(myseqs,routeU->vehicle);
			if (!y->estUnDepot) // exchange U and (V,Y or Y,V)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq12);
				addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place);
				myseqs.push_back(ysuiv->seqi_n);
				resultMoves[0][2] = seq->evaluationLB(myseqs,routeU->vehicle);
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq21);
				addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place);
				myseqs.push_back(ysuiv->seqi_n);
				resultMoves[0][3] = seq->evaluationLB(myseqs,routeU->vehicle);
			}
		}

		if (!turned)
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			addSeqDataInPieces(x,noeudV->place-1-x->place);
			myseqs.push_back(noeudU->seq1);
			myseqs.push_back(noeudV->seqi_n);
			resultMoves[1][0] = seq->evaluationLB(myseqs,routeU->vehicle);
		}
		if (!noeudV->estUnDepot) // exchange U and V
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			myseqs.push_back(noeudV->seq1);
			addSeqDataInPieces(x,noeudV->place-1-x->place);
			myseqs.push_back(noeudU->seq1);
			myseqs.push_back(y->seqi_n);
			resultMoves[1][1] = seq->evaluationLB(myseqs,routeU->vehicle);
			if (!y->estUnDepot && turned) // exchange U and (V,Y or Y,V)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq12);
				addSeqDataInPieces(x,noeudV->place-1-x->place);
				myseqs.push_back(noeudU->seq1);
				myseqs.push_back(ysuiv->seqi_n);
				resultMoves[1][2] = seq->evaluationLB(myseqs,routeU->vehicle);

				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq21);
				addSeqDataInPieces(x,noeudV->place-1-x->place);
				myseqs.push_back(noeudU->seq1);
				myseqs.push_back(ysuiv->seqi_n);
				resultMoves[1][3] = seq->evaluationLB(myseqs,routeU->vehicle);
			}
		}

		if (!x->estUnDepot)
		{
			if (!turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
				myseqs.push_back(noeudU->seq12);
				myseqs.push_back(noeudV->seqi_n);
				resultMoves[2][0] = seq->evaluationLB(myseqs,routeU->vehicle);
			}

			if (!noeudV->estUnDepot) // exchange U,X and V
			{
				if (!turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq1);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
					myseqs.push_back(noeudU->seq12);
					myseqs.push_back(y->seqi_n);
					resultMoves[2][1] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
				if (!y->estUnDepot) // exchange U,X and (V,Y or Y,V)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq12);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
					myseqs.push_back(noeudU->seq12);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[2][2] = seq->evaluationLB(myseqs,routeU->vehicle);

					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq21);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
					myseqs.push_back(noeudU->seq12);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[2][3] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
			}

			if (!turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
				myseqs.push_back(noeudU->seq21);
				myseqs.push_back(noeudV->seqi_n);
				resultMoves[3][0] = seq->evaluationLB(myseqs,routeU->vehicle);
			}
			if (!noeudV->estUnDepot) // exchange X,U and V
			{
				if (!turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq1);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
					myseqs.push_back(noeudU->seq21);
					myseqs.push_back(y->seqi_n);
					resultMoves[3][1] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
				if (!y->estUnDepot) // exchange X,U and (V,Y or Y,V)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq12);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
					myseqs.push_back(noeudU->seq21);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[3][2] = seq->evaluationLB(myseqs,routeU->vehicle);

					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq21);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
					myseqs.push_back(noeudU->seq21);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[3][3] = seq->evaluationLB(myseqs,routeU->vehicle);
				}
			}
		}
	}
	else if (decalage == 2) // U and V are almost consecutive, we are in this situation : UXVY. Useful moves must be tested as special cases
	{
		// XUVY
		resultMoves[1][0] = seq->evaluationLB(Upred->seq0_i,noeudU->seq21,noeudV->seqi_n,routeU->vehicle);
		if (!noeudV->estUnDepot)
		{
			// VXUY
			resultMoves[1][1] = seq->evaluationLB(Upred->seq0_i,noeudV->seq1,noeudU->seq21,y->seqi_n,routeU->vehicle);
			// VUXY
			resultMoves[2][1] = seq->evaluationLB(Upred->seq0_i,noeudV->seq1,noeudU->seq12,y->seqi_n,routeU->vehicle);
			if (!y->estUnDepot)
			{
				// VYXU
				resultMoves[1][2] = seq->evaluationLB(Upred->seq0_i,noeudV->seq12,noeudU->seq21,ysuiv->seqi_n,routeU->vehicle);
				// YVXU
				resultMoves[1][3] = seq->evaluationLB(Upred->seq0_i,noeudV->seq21,noeudU->seq21,ysuiv->seqi_n,routeU->vehicle);
				// VYUX
				resultMoves[2][2] = seq->evaluationLB(Upred->seq0_i,noeudV->seq12,noeudU->seq12,ysuiv->seqi_n,routeU->vehicle);
				// YVUX
				resultMoves[2][3] = seq->evaluationLB(Upred->seq0_i,noeudV->seq21,noeudU->seq12,ysuiv->seqi_n,routeU->vehicle);
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
		if (!noeudV->estUnDepot && turned)
		{

			if (shouldBeTested[0][1])
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq1);
				addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place);
				myseqs.push_back(y->seqi_n);
				resultMoves[0][1] = seq->evaluation(myseqs,routeU->vehicle);
			}

			if (!y->estUnDepot) // exchange U and (V,Y or Y,V)
			{
				if (shouldBeTested[0][2])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq12);
					addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[0][2] = seq->evaluation(myseqs,routeU->vehicle);
				}

				if (shouldBeTested[0][3])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq21);
					addSeqDataInPieces(noeudU,noeudV->place-1-noeudU->place);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[0][3] = seq->evaluation(myseqs,routeU->vehicle);
				}
			}
		}

		if (!turned && shouldBeTested[1][0])
		{
			myseqs.clear();
			myseqs.push_back(Upred->seq0_i);
			addSeqDataInPieces(x,noeudV->place-1-x->place);
			myseqs.push_back(noeudU->seq1);
			myseqs.push_back(noeudV->seqi_n);
			resultMoves[1][0] = seq->evaluation(myseqs,routeU->vehicle);
		}
		if (!noeudV->estUnDepot) // exchange U and V
		{
			if (shouldBeTested[1][1])
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				myseqs.push_back(noeudV->seq1);
				addSeqDataInPieces(x,noeudV->place-1-x->place);
				myseqs.push_back(noeudU->seq1);
				myseqs.push_back(y->seqi_n);
				resultMoves[1][1] = seq->evaluation(myseqs,routeU->vehicle);
			}

			if (!y->estUnDepot && turned) // exchange U and (V,Y or Y,V)
			{
				if (shouldBeTested[1][2])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq12);
					addSeqDataInPieces(x,noeudV->place-1-x->place);
					myseqs.push_back(noeudU->seq1);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[1][2] = seq->evaluation(myseqs,routeU->vehicle);
				}

				if (shouldBeTested[1][3])
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq21);
					addSeqDataInPieces(x,noeudV->place-1-x->place);
					myseqs.push_back(noeudU->seq1);
					myseqs.push_back(ysuiv->seqi_n);
					resultMoves[1][3] = seq->evaluation(myseqs,routeU->vehicle);
				}
			}
		}

		if (!x->estUnDepot)
		{
			if (shouldBeTested[2][0] && !turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
				myseqs.push_back(noeudU->seq12);
				myseqs.push_back(noeudV->seqi_n);
				resultMoves[2][0] = seq->evaluation(myseqs,routeU->vehicle);
			}

			if (!noeudV->estUnDepot) // exchange U,X and V
			{
				if (shouldBeTested[2][1] && !turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq1);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
					myseqs.push_back(noeudU->seq12);
					myseqs.push_back(y->seqi_n);
					resultMoves[2][1] = seq->evaluation(myseqs,routeU->vehicle);
				}
				if (!y->estUnDepot) // exchange U,X and (V,Y or Y,V)
				{
					if (shouldBeTested[2][2])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(noeudV->seq12);
						addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
						myseqs.push_back(noeudU->seq12);
						myseqs.push_back(ysuiv->seqi_n);
						resultMoves[2][2] = seq->evaluation(myseqs,routeU->vehicle);
					}

					if (shouldBeTested[2][3])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(noeudV->seq21);
						addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
						myseqs.push_back(noeudU->seq12);
						myseqs.push_back(ysuiv->seqi_n);
						resultMoves[2][3] = seq->evaluation(myseqs,routeU->vehicle);
					}
				}
			}

			if (shouldBeTested[3][0] && !turned)
			{
				myseqs.clear();
				myseqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
				myseqs.push_back(noeudU->seq21);
				myseqs.push_back(noeudV->seqi_n);
				resultMoves[3][0] = seq->evaluation(myseqs,routeU->vehicle);
			}
			if (!noeudV->estUnDepot) // exchange X,U and V
			{
				if (shouldBeTested[3][1] && !turned)
				{
					myseqs.clear();
					myseqs.push_back(Upred->seq0_i);
					myseqs.push_back(noeudV->seq1);
					addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
					myseqs.push_back(noeudU->seq21);
					myseqs.push_back(y->seqi_n);
					resultMoves[3][1] = seq->evaluation(myseqs,routeU->vehicle);
				}
				if (!y->estUnDepot) // exchange X,U and (V,Y or Y,V)
				{
					if (shouldBeTested[3][2])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(noeudV->seq12);
						addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
						myseqs.push_back(noeudU->seq21);
						myseqs.push_back(ysuiv->seqi_n);
						resultMoves[3][2] = seq->evaluation(myseqs,routeU->vehicle);
					}

					if (shouldBeTested[3][3])
					{
						myseqs.clear();
						myseqs.push_back(Upred->seq0_i);
						myseqs.push_back(noeudV->seq21);
						addSeqDataInPieces(xsuiv,noeudV->place-1-xsuiv->place);
						myseqs.push_back(noeudU->seq21);
						myseqs.push_back(ysuiv->seqi_n);
						resultMoves[3][3] = seq->evaluation(myseqs,routeU->vehicle);
					}
				}
			}
		}
	}
	else if (decalage == 2)
	{
		// XUVY
		if (shouldBeTested[1][0]) resultMoves[1][0] = seq->evaluation(Upred->seq0_i,noeudU->seq21,noeudV->seqi_n,routeU->vehicle);
		if (!noeudV->estUnDepot)
		{
			// VXUY
			if (shouldBeTested[1][1]) resultMoves[1][1] = seq->evaluation(Upred->seq0_i,noeudV->seq1,noeudU->seq21,y->seqi_n,routeU->vehicle);
			// VUXY
			if (shouldBeTested[2][1]) resultMoves[2][1] = seq->evaluation(Upred->seq0_i,noeudV->seq1,noeudU->seq12,y->seqi_n,routeU->vehicle);
			if (!y->estUnDepot)
			{
				// VYXU
				if (shouldBeTested[1][2]) resultMoves[1][2] = seq->evaluation(Upred->seq0_i,noeudV->seq12,noeudU->seq21,ysuiv->seqi_n,routeU->vehicle);
				// YVXU
				if (shouldBeTested[1][3]) resultMoves[1][3] = seq->evaluation(Upred->seq0_i,noeudV->seq21,noeudU->seq21,ysuiv->seqi_n,routeU->vehicle);
				// VYUX
				if (shouldBeTested[2][2]) resultMoves[2][2] = seq->evaluation(Upred->seq0_i,noeudV->seq12,noeudU->seq12,ysuiv->seqi_n,routeU->vehicle);
				// YVUX
				if (shouldBeTested[2][3]) resultMoves[2][3] = seq->evaluation(Upred->seq0_i,noeudV->seq21,noeudU->seq12,ysuiv->seqi_n,routeU->vehicle);
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
		noeudU = tempU ;
		noeudV = tempV ;
		x = noeudU->suiv ;
		y = noeudV->suiv ;
		return 0 ;
	}


	///////////////////////////////////////////
	// HERE STARTS THE MOVE APPLICATION PROCESS
	///////////////////////////////////////////

	reinitMoves(routeU);  // Say that all moves involving this route must be tested again

	// Update the solution
	if (decalage >= 3)
	{
		if ( ibest == 1 || ibest == 2) insertNoeud(noeudU,placeV);
		if ( ibest == 2 ) insertNoeud(x,noeudU);
		if ( ibest == 3 ) { insertNoeud(x,placeV); insertNoeud(noeudU,x); }

		if ( jbest == 1 || jbest == 2) insertNoeud(noeudV,placeU);
		if ( jbest == 2) insertNoeud(y,noeudV);
		if ( jbest == 3 ) { insertNoeud(y,placeU); insertNoeud(noeudV,y); }
	}

	// Special cases of decalage == 2
	else if (decalage == 2)
	{
		if (ibest == 1 && jbest == 0) { insertNoeud(noeudU,x); }
		else if (ibest == 1 && jbest == 1) { insertNoeud(x,noeudV); insertNoeud(noeudU,x);  }
		else if (ibest == 1 && jbest == 2) { insertNoeud(x,y); insertNoeud(noeudU,x);  }
		else if (ibest == 1 && jbest == 3) { insertNoeud(noeudV,y); insertNoeud(x,noeudV); insertNoeud(noeudU,x);  }
		else if (ibest == 2 && jbest == 1) { insertNoeud(noeudV,noeudU->pred) ;  }
		else if (ibest == 2 && jbest == 2) { insertNoeud(noeudV,noeudU->pred) ;  insertNoeud(y,noeudV) ; }
		else if (ibest == 2 && jbest == 3) { insertNoeud(y,noeudU->pred) ;  insertNoeud(noeudV,y) ; }
		else throw string ("ERROR move intra-route") ;
	} 

	routeU->updateRouteData(false); // Update the pre-processed data on the subsequences of the route
	setRouteVide(); // Keep a pointer on the first empty route
	rechercheTerminee = false ; // Not finished the search
	nbIntraSwap ++ ;
	noeudU = tempU ;
	noeudV = tempV ;
	x = noeudU->suiv ;
	y = noeudV->suiv ;
	return 1 ; // Return Success
}


int LocalSearch::intraRoute2Opt ()
{
	// Evaluation procedure for 2-Opt
	Noeud * nodeNum = noeudU->suiv ;
	Noeud * nodeUpred = noeudU->pred ;
	Noeud * temp ;
	SeqData * seq = noeudU->seq0_i ;

	double cost ;
	double costZero = routeU->currentRouteCost ;

	myseqs.clear();
	myseqs.push_back(noeudU->pred->seq0_i);
	addReverseSeqDataInPieces(noeudU,noeudV->place-noeudU->place);
	myseqs.push_back(noeudV->suiv->seqi_n);

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
	noeudU->pred = nodeNum ;
	noeudU->suiv = y ;

	while ( nodeNum != noeudV )
	{
		temp = nodeNum->suiv ;
		nodeNum->suiv = nodeNum->pred ;
		nodeNum->pred = temp ;
		nodeNum = temp ;
	}

	noeudV->suiv = noeudV->pred ;
	noeudV->pred = nodeUpred ;
	nodeUpred->suiv = noeudV ;
	y->pred = noeudU ;
	routeU->updateRouteData(false);  // Update the pre-processed data on the subsequences of the route
	setRouteVide(); // Keep a pointer on the first empty route
	rechercheTerminee = false ; // Not finished the search
	nbIntra2Opt ++ ;
	return 1 ; // Return Success
}

void LocalSearch::melangeParcours ()
{
	// Shuffling the ordreParcours vector
	int j,temp ;
  for (int i = 0 ; i < (int)ordreParcours.size() - 1 ; i++)
  {
    j = i + rand() % ((int)ordreParcours.size() - i) ;
    temp = ordreParcours[i] ;
    ordreParcours[i] = ordreParcours[j] ;
    ordreParcours[j] = temp ;
  }
}

void LocalSearch::updateMoves ()
{
	int client, client2, size ;
  for (int i=0 ; i<( int)ordreParcours.size() ; i++)
  {
    client = ordreParcours[i] ;
    clients[client].moves.clear();
    size = (int)params->cli[client].inverseCloseCustomers.size();
    for (int a1 = 0 ; a1 < size ; a1++ )
    {
      client2 = params->cli[client].inverseCloseCustomers[a1] ;
      if (client2 >= params->nbDepots && clients[client2].estPresent) clients[client].moves.push_back(client2);
    }
  }
}

void LocalSearch::setRouteVide()
{
	routeVide = NULL ;
	int route = 0 ;
	while (routeVide == NULL && route < params->nbVehiculesPerDep )
	{
		if (depots[route].suiv->estUnDepot) routeVide = depots[route].route ;
		route ++ ;
	}
}

void LocalSearch::nodeTestedForEachRoute (int cli)
{
	for (int route = 0 ; route < (int)params->nombreVehicules ; route ++)
		routes[route].nodeAndRouteTested[cli]=true ;	
}

void LocalSearch::insertNoeud(Noeud * U, Noeud * V)
{
	if (U->pred != V && U != V)
	{
		U->pred->suiv = U->suiv ;
		U->suiv->pred = U->pred ;
		V->suiv->pred = U ;
		U->pred = V ;
		U->suiv = V->suiv ;
		V->suiv = U ;
		U->route = V->route ;
	}
}

void LocalSearch::removeNoeud(Noeud * U)
{
	Noeud * temp =  U->suiv ;
	reinitMoves(U->route);
	U->pred->suiv = U->suiv ;
	U->suiv->pred = U->pred ;
	U->route = NULL ;
	U->pred = NULL ;
	U->suiv = NULL ;
	temp->route->updateRouteData(false);
	setRouteVide();

	// Managing the other data structures
	individu->chromT.pop_back();
	removeOP(U->customer_or_depot_idx);
	U->estPresent = false ;
}

void LocalSearch::addNoeud(Noeud * U, Noeud * V)
{
	reinitMoves(V->route);
	V->suiv->pred = U ;
	U->pred = V;
	U->suiv = V->suiv ;
	V->suiv = U ;

	// Update the routes
	U->estPresent = true ;
	U->route = V->route ;
	U->route->updateRouteData(false);
	setRouteVide();

	// Manage the other data structures
	individu->chromT.push_back(0);
	addOP(U->customer_or_depot_idx);
}

void LocalSearch::removeOP (int client) 
{
	// Remove one occurence in the ordreParcours structure
	int it = 0  ;
	while (ordreParcours[it] != client) { it ++ ; }
	ordreParcours[it] = ordreParcours[(int)ordreParcours.size() - 1];
	ordreParcours.pop_back();
}

void LocalSearch::addOP (int client) 
{
	// Add one element in the ordreParcours structure
	int it, temp2 ;
	if (ordreParcours.size() != 0)
	{
		it = (int)rand() % ordreParcours.size() ;
		temp2 = ordreParcours[it] ;
		ordreParcours[it] = client ;
		ordreParcours.push_back(temp2);
	}
	else
		ordreParcours.push_back(client);
}

void LocalSearch::reinitMoves(Route * r)
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
		r->nodeAndRouteTested[i] = false ;

	for (Noeud * tempNoeud = r->depot->suiv ; !tempNoeud->estUnDepot ; tempNoeud = tempNoeud->suiv)
		for (int route = 0 ; route < params->nbVehiculesPerDep ; route ++)
			depots[route].route->nodeAndRouteTested[tempNoeud->customer_or_depot_idx] = false ;
}

void LocalSearch::reinitAllMoves()
{
  for (int route = 0 ; route < params->nbVehiculesPerDep ; route ++)
    for (int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
      depots[route].route->nodeAndRouteTested[i] = false ;
}

void LocalSearch::addSeqDataInPieces (Noeud * node, int length)
{
	Noeud * curNode = node ;
	SeqData * curSeq ;
	int curLength = length + 1 ;
	int temp ;
	while (curLength > 0)
	{
		temp = min(curLength,params->sizeSD) ;
		curSeq = curNode->seqi_j[temp-1];
		myseqs.push_back(curSeq);
		curLength -= temp ;
		if (curLength != 0) curNode = clients[curSeq->lastNode].suiv ;
	}
}

void LocalSearch::addReverseSeqDataInPieces (Noeud * node, int length)
{
	Noeud * curNode = node ;
	SeqData * curSeq ;
	list <SeqData*> myseqTemp ;
	int curLength = length + 1 ;
	int temp ;
	while (curLength > 0)
	{
		temp = min(curLength,params->sizeSD) ;
		curSeq = curNode->seqj_i[temp-1];
		myseqTemp.push_back(curSeq);
		curLength -= temp ;
		if (curLength != 0) curNode = clients[curSeq->firstNode].suiv ;
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
	vector < Noeud * > tempNoeud ;
	vector <int> temp2 ;

	nbVeh = params->nombreVehicules ;
	clients = new Noeud [params->nbClients + params->nbDepots + 1] ;
	depots = new Noeud [nbVeh] ;
	depotsFin = new Noeud [nbVeh] ;
	routes = new Route [nbVeh] ; 

	for (int i = 0 ; i <  params->nbClients + params->nbDepots ; i ++ )
		clients[i] = Noeud(false,i,false,NULL,NULL,NULL,params);

	for (int i = 0 ; i < nbVeh ; i++ )
	{
		depots[i] = Noeud(true,params->ordreVehicules[i].depotNumber,false,NULL,NULL,NULL,params);
		depotsFin[i] = Noeud(true,params->ordreVehicules[i].depotNumber,false,NULL,NULL,NULL,params);
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

  routeVide = NULL;

	for (int i=params->nbDepots ; i < params->nbDepots + params->nbClients ; i++)
		ordreParcours.push_back(i);

	// Initialization of the SeqData structures, for all nodes
	// These structures were designed to stay contiguous in memory (however it had only a negligible impact on performance)
	int nbSeqsSet = 0 ;
	int taillemyseqDatas = (params->sizeSD + params->sizeSD + 4)*(params->nbClients + params->nbDepots + 2*params->nbVehiculesPerDep) ;
	SeqData * myseqDatas = new SeqData [taillemyseqDatas] ;
	seqdeb = myseqDatas ;
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i++) 
	{
		for (int ii=nbSeqsSet ; ii < nbSeqsSet+4+params->sizeSD+params->sizeSD ; ii++)
			myseqDatas[ii].initialisation(clients[i].customer_or_depot_idx,params,individu,false);  

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
			myseqDatas[ii].initialisation(depots[i].customer_or_depot_idx,params,individu,false);

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
			myseqDatas[ii].initialisation(depotsFin[i].customer_or_depot_idx,params,individu,false);

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

