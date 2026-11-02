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

#ifndef SEQDATA_H
#define SEQDATA_H 

#include "Params.h" 
class Individu ;

class SeqData
{
public:

	Params * params ; // Access to the problem parameters
	int firstNode ; // first node of the SeqData
	int lastNode ; // last node of the SeqData
	double load ; // total demand on the SeqData
	double cost ; // total cost of the SeqData
	double distance ; // total distance of the SeqData

	// bestCostNodes[i] contains the sequence of nodes when starting with the depot 0 and finishing with node i
	vector < vector < vector <int> > > bestCostNodes ;

	// More general version, which can deal with an unlimited number of service modes
	// Used for problems with turn penalties
	// bestCost[i][j] gives the best cost 
	// when starting the first service with its mode i
	// and finishing the last service with its mode j
	vector < vector < double > > bestCost ;

	// two temporary vectors, initialized with the good size, and used to perform evaluations
	vector < double > distanceTemp ;
	vector < double > distanceTemp2 ;

	// just to say if the structures are already created or not.
	bool isInitialized ;

	// Construction operators for data pre-processing
	void initialisation(int U_customer_or_depot_idx, Params * mesParams, Individu * myIndiv, bool isForPathTracking);
	void concatOneAfter(SeqData * seq,int V_customer_or_depot_idx, Individu * myIndiv);
	void concatOneAfterWithPathTracking(SeqData * seq,int V_customer_or_depot_idx, Individu * myIndiv); // used to track the path when printing the final solution
	void concatOneBefore(SeqData * seq,int V_customer_or_depot_idx, Individu * myIndiv);
	
	// Route evaluation evaluators
	double evaluation(SeqData * seq1, Vehicle * vehicle);
	double evaluation(SeqData * seq1, SeqData * seq2, Vehicle * vehicle);
	double evaluation(SeqData * seq1, SeqData * seq2, Vehicle * vehicle, double & mydist, double & mytminex, double & myloadex);
	double evaluation(SeqData * seq1, SeqData * seq2, SeqData * seq3, Vehicle * vehicle);
	double evaluation(SeqData * seq1, SeqData * seq2, SeqData * seq3, SeqData * seq4, Vehicle * vehicle);
	double evaluation(vector <SeqData *> seqs, Vehicle * vehicle);

	// The same evaluators, but to get lower bounds on move evaluations
	double evaluationLB(SeqData * seq1, Vehicle * vehicle);
	double evaluationLB(SeqData * seq1, SeqData * seq2, Vehicle * vehicle);
	double evaluationLB(SeqData * seq1, SeqData * seq2, SeqData * seq3, Vehicle * vehicle);
	double evaluationLB(SeqData * seq1, SeqData * seq2, SeqData * seq3, SeqData * seq4, Vehicle * vehicle);
	double evaluationLB(vector <SeqData *> seqs, Vehicle * vehicle);

	// Constructors
	SeqData(Params * params);
	SeqData();

	// Destructors
	~SeqData();
};

#endif
