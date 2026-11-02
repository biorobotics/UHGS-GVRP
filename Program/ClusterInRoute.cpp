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

#include "ClusterInRoute.h"

ClusterInRoute::ClusterInRoute(void) {}

ClusterInRoute::ClusterInRoute(bool isDepot, int cluster_idx,bool isPresent, ClusterInRoute * successor , ClusterInRoute * pred, Route * route,Params * params) 
: params(params), isDepot(isDepot), cluster_idx(cluster_idx), isPresent(isPresent), successor(successor), pred(pred), route(route)
{
	place = -1 ;
}

ClusterInRoute::ClusterInRoute(ClusterInRoute const& copy)
{
	// Copy constructor
	isDepot = copy.isDepot ;
	cluster_idx = copy.cluster_idx ;
	place = copy.place ;
	isPresent = copy.isPresent ;
	successor = copy.successor ;
	pred = copy.pred ;
	route = copy.route ;
	params = copy.params ;
	moves = copy.moves ;
}

ClusterInRoute& ClusterInRoute::operator=(ClusterInRoute const& copy)
{
	// Copy constructor
	isDepot = copy.isDepot ;
	cluster_idx = copy.cluster_idx ;
	place = copy.place ;
	isPresent = copy.isPresent ;
	successor = copy.successor ;
	pred = copy.pred ;
	route = copy.route ;
	params = copy.params ;
	moves = copy.moves ;
	return *this;
}

ClusterInRoute::~ClusterInRoute(void){}

void ClusterInRoute::setRemaining()
{
	// seq1 has exactly the same meaning than seqi_j[0], but its more convenient to use and read
	seq1 = seqi_j[0] ;
	seq12 = seqi_j[1] ;
	seq123 = seqi_j[2] ;
	seq21 = seqj_i[1] ;
	seq321 = seqj_i[2] ;
}
