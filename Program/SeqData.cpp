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

#include "SeqData.h"
#include "Individu.h" 

/* --------------------------------------------------------- */
/* ROUTE EVALUATION PROCEDURES FOR GVRP */
/* --------------------------------------------------------- */

SeqData::SeqData(Params * params)
{
	this->params = params ;
	firstNode = -1 ;
	lastNode = -1 ;

	distanceTemp = vector < double > (params->ar_maxNbNodesPerClient);
	distanceTemp2 = vector < double > (params->ar_maxNbNodesPerClient);

	bestCost = vector < vector < double > > (params->ar_maxNbNodesPerClient) ;
	for (int i=0 ; i < params->ar_maxNbNodesPerClient ; i++)
		bestCost[i] =  vector < double > (params->ar_maxNbNodesPerClient);

	bestCostNodes = vector < vector < vector <int> > > (params->ar_maxNbNodesPerClient) ;
	for (int i=0 ; i < params->ar_maxNbNodesPerClient ; i++)
		bestCostNodes[i] = vector < vector <int> > (params->ar_maxNbNodesPerClient);

	isInitialized = true ;
}

SeqData::SeqData()
{
	firstNode = -1 ;
	lastNode = -1 ;
	isInitialized = false ;
}

SeqData::~SeqData()
{}

void SeqData::initialisation(int Ucour, Params * mesParams, Individu * myIndiv, bool isForPathTracking)
{
	params = mesParams ;
	load = params->cli[Ucour].demand;

	if (!isInitialized)
	{
		bestCost = vector < vector <double> > (params->ar_maxNbNodesPerClient) ;
		distanceTemp = vector < double > (params->ar_maxNbNodesPerClient);
		distanceTemp2 = vector < double > (params->ar_maxNbNodesPerClient);
		for (int i=0 ; i < params->ar_maxNbNodesPerClient ; i++)
			bestCost[i] =  vector <double> (params->ar_maxNbNodesPerClient);
		bestCostNodes = vector < vector < vector <int> > > (params->ar_maxNbNodesPerClient) ;
		for (int i=0 ; i < params->ar_maxNbNodesPerClient ; i++)
			bestCostNodes[i] = vector < vector <int> > (params->ar_maxNbNodesPerClient);
		isInitialized = true ;
	}

	for (int i=0 ; i < params->cli[Ucour].ar_nbNodes ; i++)
	{
		for (int j=0 ; j < params->cli[Ucour].ar_nbNodes ; j ++)
		{
			if (i == j) 
				bestCost[i][j] = 0.0 ;
			else 
				bestCost[i][j] = 1.e30 ;
		}
	}

	if (isForPathTracking)
	{
		for (int i=0 ; i < params->cli[Ucour].ar_nbNodes ; i++)
		{
			for (int j=0 ; j < params->cli[Ucour].ar_nbNodes ; j ++)
			{
				bestCostNodes[i][j].clear();
				if (i == j) 
					bestCostNodes[i][j].push_back(params->cli[Ucour].ar_Nodes[i]);
			}
		}
	}

	distance = 0 ;
	firstNode = Ucour ;
	lastNode = Ucour ;
}

void SeqData::concatOneAfter(SeqData * seq, int Vcour, Individu * myIndiv)
{
	double tempc ;
	Client * firstCli = &params->cli[seq->firstNode] ;
	Client * lastCli  = &params->cli[seq->lastNode] ;
	Client * vCourCli = &params->cli[Vcour] ;
	distance = 1.e20 ;

	for (int i=0 ; i < firstCli->ar_nbNodes ; i++)
	{
		for (int j=0 ; j < vCourCli->ar_nbNodes ; j++)
		{
			bestCost[i][j] = 1.e30 ;
			for (int k=0 ; k < lastCli->ar_nbNodes ; k++)
			{
				tempc = seq->bestCost[i][k] + params->cost_mat[lastCli->ar_Nodes[k]][vCourCli->ar_Nodes[j]] ;
				if (tempc < bestCost[i][j]) bestCost[i][j] = tempc;
				if (tempc < distance) distance = tempc ;
			}
		}
	}

	load = seq->load + params->cli[Vcour].demand;
	firstNode = seq->firstNode ;
	lastNode = Vcour ;
}

void SeqData::concatOneAfterWithPathTracking(SeqData * seq, int Vcour, Individu * myIndiv)
{
	double tempc ;
	Client * firstCli = &params->cli[seq->firstNode] ;
	Client * lastCli  = &params->cli[seq->lastNode] ;
	Client * vCourCli = &params->cli[Vcour] ;
	distance = 1.e20 ;

	for (int i=0 ; i < firstCli->ar_nbNodes ; i++)
	{
		for (int j=0 ; j < vCourCli->ar_nbNodes ; j++)
		{
			bestCost[i][j] = 1.e30 ;
			int myK = -1;
			for (int k=0 ; k < lastCli->ar_nbNodes ; k++)
			{
				tempc = seq->bestCost[i][k] + params->cost_mat[lastCli->ar_Nodes[k]][vCourCli->ar_Nodes[j]] ;
				if (tempc < bestCost[i][j]) 
				{
					bestCost[i][j] = tempc ;
					myK = k ;
				}
				if (tempc < distance) distance = tempc ;
			}

			if (myK == -1)
			{
				cout << "Issue undetermined K during Label propagation" << endl ;
				throw string ("Issue undetermined K during Label propagation");
			}
			bestCostNodes[i][j].clear();
			bestCostNodes[i][j] = seq->bestCostNodes[i][myK] ;
			bestCostNodes[i][j].push_back(params->cli[Vcour].ar_Nodes[j]);
		}
	}

	load = seq->load + params->cli[Vcour].demand;
	firstNode = seq->firstNode ;
	lastNode = Vcour ;
}


void SeqData::concatOneBefore(SeqData * seq,int Vcour, Individu * myIndiv)
{ 
	double tempc ;
	Client * firstCli = &params->cli[seq->firstNode] ;
	Client * lastCli  = &params->cli[seq->lastNode] ;
	Client * vCourCli = &params->cli[Vcour] ;
	distance = 1.e20 ;

	for (int i=0 ; i < vCourCli->ar_nbNodes ; i++)
	{
		for (int j=0 ; j < lastCli->ar_nbNodes ; j++)
		{
			bestCost[i][j] = 1.e30 ;
			for (int k=0 ; k < firstCli->ar_nbNodes ; k++)
			{
				tempc = params->cost_mat[vCourCli->ar_Nodes[i]][firstCli->ar_Nodes[k]] + seq->bestCost[k][j] ;
				if (tempc < bestCost[i][j]) bestCost[i][j] = tempc ;
				if (tempc < distance) distance = tempc ;
			}
		}
	}

	load = seq->load + params->cli[Vcour].demand;
	firstNode = Vcour ;
	lastNode = seq->lastNode ;
}

// Returns cost of seq1
double SeqData::evaluation(SeqData * seq1, Vehicle * vehicle) 
{
	return seq1->bestCost[0][0] + max(seq1->load - vehicle->vehicleCapacity,0.0)*params->penalityCapa ;
}

// Computes cost of concatenating seq1 to seq2 (seq2 comes after seq1)
double SeqData::evaluation(SeqData * seq1, SeqData * seq2, Vehicle * vehicle) 
{
	Client * cli1 = &params->cli[seq1->lastNode] ;
	Client * cli2 = &params->cli[seq2->firstNode] ;
	double bestc = 1.e30 ;
	double tempc ;

	for (int i=0 ; i < cli1->ar_nbNodes ; i++)
	{
		for (int j=0 ; j < cli2->ar_nbNodes ; j++)
		{
			tempc = seq1->bestCost[0][i] +
				params->cost_mat[cli1->ar_Nodes[i]][cli2->ar_Nodes[j]]
			+ seq2->bestCost[j][0] ;
			if (tempc < bestc) bestc = tempc ;
		}
	}

	return bestc + max(seq1->load + seq2->load - vehicle->vehicleCapacity,0.0)*params->penalityCapa ;
}

double SeqData::evaluation(SeqData * seq1, SeqData * seq2, Vehicle * vehicle, double & mydist, double & mytminex, double & myloadex) 
{
	Client * cli1 = &params->cli[seq1->lastNode] ;
	Client * cli2 = &params->cli[seq2->firstNode] ;
	double bestc = 1.e30 ;
	double tempc ;

	for (int i=0 ; i < cli1->ar_nbNodes ; i++)
	{
		for (int j=0 ; j < cli2->ar_nbNodes ; j++)
		{
			tempc = seq1->bestCost[0][i] +
				params->cost_mat[cli1->ar_Nodes[i]][cli2->ar_Nodes[j]]
			+ seq2->bestCost[j][0] ;
			if (tempc < bestc) bestc = tempc ;
		}
	}

	mydist = bestc ;
	mytminex = 0 ;
	myloadex = max(seq1->load + seq2->load - vehicle->vehicleCapacity,0.0);

	return mydist + myloadex*params->penalityCapa ;
}

// Compute cost of concatenating seq1, seq2, seq3 in that order
double SeqData::evaluation(SeqData * seq1, SeqData * seq2, SeqData * seq3, Vehicle * vehicle) 
{
	Client * cliLast1 = &params->cli[seq1->lastNode] ;
	Client * cliLast2 = &params->cli[seq2->lastNode] ;
	Client * cliFirst2 = &params->cli[seq2->firstNode] ;
	Client * cliFirst3 = &params->cli[seq3->firstNode] ;

	double bestc ;
	double tempc ;

  // Compute cost-to-come to each mode at cliFirst2
	for (int i=0 ; i < cliFirst2->ar_nbNodes ; i++)
	{
		distanceTemp2[i] = 1.e30 ;
		for (int k=0 ; k < cliLast1->ar_nbNodes ; k++)
		{
			tempc = seq1->bestCost[0][k] + params->cost_mat[cliLast1->ar_Nodes[k]][cliFirst2->ar_Nodes[i]] ;
			if (tempc < distanceTemp2[i]) distanceTemp2[i] = tempc ;
		}
	}

  // Compute cost-to-come to each mode at cliLast2
	for (int i=0 ; i < cliLast2->ar_nbNodes ; i++)
	{
		distanceTemp[i] = 1.e30 ;
		for (int k=0 ; k < cliFirst2->ar_nbNodes ; k++)
		{
			tempc = distanceTemp2[k] + seq2->bestCost[k][i] ;
			if (tempc < distanceTemp[i]) distanceTemp[i] = tempc ;
		}
	}

  // Compute cost
	bestc = 1.e30 ;
	for (int i=0 ; i < cliLast2->ar_nbNodes ; i++)
	{
		for (int j=0 ; j < cliFirst3->ar_nbNodes ; j++)
		{
			tempc = distanceTemp[i] 
			+ params->cost_mat[cliLast2->ar_Nodes[i]][cliFirst3->ar_Nodes[j]]
			+ seq3->bestCost[j][0] ;
			if (tempc < bestc) bestc = tempc ;
		}
	}

	return bestc + max(seq1->load + seq2->load  + seq3->load - vehicle->vehicleCapacity,0.0)*params->penalityCapa ;
}

// Compute cost of concatenating seq1, seq2, seq3, seq4 in that order
double SeqData::evaluation(SeqData * seq1, SeqData * seq2, SeqData * seq3, SeqData * seq4, Vehicle * vehicle) 
{
	vector <SeqData *> seqs ;
	seqs.push_back(seq1);
	seqs.push_back(seq2);
	seqs.push_back(seq3);
	seqs.push_back(seq4);
	return evaluation(seqs,vehicle);
}

// Compute cost of concatenating seqs in order
double SeqData::evaluation(vector <SeqData *> seqs, Vehicle * vehicle) 
{
	double tempc ;
	Client * cliPredLast ;
	Client * cliFirst ;
	Client * cliLast ;

	SeqData *seqbPred = seqs[0];
	SeqData *seqb = seqs[1];
	double loadTemp = seqbPred->load ;
	distanceTemp = seqbPred->bestCost[0] ;

	for (int s=1 ; s < (int)seqs.size() ; s++)
	{
		seqbPred = seqs[s-1];
		seqb = seqs[s];
		cliPredLast = &params->cli[seqbPred->lastNode] ;
		cliFirst = &params->cli[seqb->firstNode] ;
		cliLast = &params->cli[seqb->lastNode] ;
		loadTemp += seqb->load ;

		for (int i=0 ; i < cliFirst->ar_nbNodes ; i++)
		{
			distanceTemp2[i] = 1.e30 ;
			for (int k=0 ; k < cliPredLast->ar_nbNodes ; k++)
			{
				tempc = distanceTemp[k] + params->cost_mat[cliPredLast->ar_Nodes[k]][cliFirst->ar_Nodes[i]] ;
				if (tempc < distanceTemp2[i]) distanceTemp2[i] = tempc ;
			}
		}

		for (int i=0 ; i < cliLast->ar_nbNodes ; i++)
		{
			distanceTemp[i] = 1.e30 ;
			for (int k=0 ; k < cliFirst->ar_nbNodes ; k++)
			{
				tempc = distanceTemp2[k] + seqb->bestCost[k][i] ;
				if (tempc < distanceTemp[i]) distanceTemp[i] = tempc ;
			}
		}		
	}
	return distanceTemp[0] + max(loadTemp - vehicle->vehicleCapacity,0.0)*params->penalityCapa ;
}

double SeqData::evaluationLB(SeqData * seq1, Vehicle * vehicle) 
{
	return seq1->distance + max(seq1->load - vehicle->vehicleCapacity,0.0)*params->penalityCapa ;
}

double SeqData::evaluationLB(SeqData * seq1, SeqData * seq2, Vehicle * vehicle) 
{
	return seq1->distance + seq2->distance + params->timeCost[seq1->lastNode][seq2->firstNode]
	+ max(seq1->load + seq2->load - vehicle->vehicleCapacity,0.0)*params->penalityCapa ;
}

double SeqData::evaluationLB(SeqData * seq1, SeqData * seq2, SeqData * seq3, Vehicle * vehicle) 
{
	return seq1->distance + seq2->distance + seq3->distance + params->timeCost[seq1->lastNode][seq2->firstNode] + params->timeCost[seq2->lastNode][seq3->firstNode]
	+ max(seq1->load + seq2->load + seq3->load - vehicle->vehicleCapacity,0.0)*params->penalityCapa ;
}

double SeqData::evaluationLB(SeqData * seq1, SeqData * seq2, SeqData * seq3, SeqData * seq4, Vehicle * vehicle) 
{
	return seq1->distance + seq2->distance + seq3->distance + seq4->distance +params->timeCost[seq1->lastNode][seq2->firstNode] 
	+ params->timeCost[seq2->lastNode][seq3->firstNode]
	+ params->timeCost[seq3->lastNode][seq4->firstNode]
	+ max(seq1->load + seq2->load + seq3->load + seq4->load - vehicle->vehicleCapacity,0.0)*params->penalityCapa ;
}

double SeqData::evaluationLB(vector <SeqData *> seqs, Vehicle * vehicle) 
{
	SeqData *seqbPred = seqs[0];
	SeqData *seqb = seqs[1];
	double myLoadTemp = seqbPred->load ;
	double myDistanceTemp = seqbPred->distance ;

	for (int i=1 ; i<(int)seqs.size() ; i++)
	{
		seqbPred = seqs[i-1];
		seqb = seqs[i];
		myDistanceTemp += params->timeCost[seqbPred->lastNode][seqb->firstNode] + seqb->distance ;
		myLoadTemp += seqb->load ;
	}
	return myDistanceTemp + max(myLoadTemp - vehicle->vehicleCapacity,0.0)*params->penalityCapa ;
}
