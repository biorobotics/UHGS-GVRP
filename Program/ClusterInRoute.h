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

#ifndef NOEUD_H
#define NOEUD_H

#include <iostream>
using namespace std;
class Route ;
#include "Route.h"

class ClusterInRoute
{

public :

// Access to the data of the problem
Params * params ;

// is-it a depot
bool isDepot ;

// index of the depot or customer
int cluster_idx ;

// place in the route
int place ;

// is this customer used
// (all customer nodes are created, but not necessarily inserted in the sequence)
bool isPresent ;

// next depot or client in the route
ClusterInRoute * successor ;

// previous depot or client in the route
ClusterInRoute * pred ;

// associated route
Route * route ;

// pointer towards the preprocessed SeqData data structures
// "i" is considered to be the current customer
vector <SeqData *> seqi_j ; // data for (i,j) with j > i
vector <SeqData *> seqj_i ; // data for (j,i) (for the same subsequence as i_j, but reversed)
SeqData * seq0_i ; // data for (0,i)
SeqData * seqi_n ; // data for (i,n), n is the end of the route
SeqData * seqi_0 ; // data for (i,0) (for the reversed route)
SeqData * seqn_i ; // data for (n,i) (for the reversed route)
// the same pointers as (i,j) for some values, but simpler to call
SeqData * seq1 ; // data for (i) 
SeqData * seq12 ; // data for (i,i+1)
SeqData * seq21 ; // data for (i+1,i)
SeqData * seq123 ; // data for (i,i+1,i+2)
SeqData * seq321 ; // data for (i+2,i+1,i)

// possible moves for this customer (granular search)
vector < int > moves ;

// constructor 1
ClusterInRoute(void);
	
// constructor 2
ClusterInRoute(bool isDepot, int customer_or_depot_index, bool isPresent, ClusterInRoute * successor , ClusterInRoute * pred, Route * route,Params * params);

// destructor
~ClusterInRoute(void);

// Copy constructor
ClusterInRoute(ClusterInRoute const& copy) ;

// Assignment operator in terms of the copy constructor
ClusterInRoute& operator=(ClusterInRoute const& copy);

// little function to correctly initialize the pointers
void setRemaining();

};

#endif
