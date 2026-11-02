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

#include "Noeud.h"

Noeud::Noeud(void) {}

Noeud::Noeud(bool isDepot, int customer_or_depot_idx,bool isPresent, Noeud * successor , Noeud * pred, Route * route,Params * params) 
: params(params), isDepot(isDepot), customer_or_depot_idx(customer_or_depot_idx), isPresent(isPresent), successor(successor), pred(pred), route(route)
{
	place = -1 ;
}

Noeud::Noeud(Noeud const& copy)
{
	// Copy constructor
	isDepot = copy.isDepot ;
	customer_or_depot_idx = copy.customer_or_depot_idx ;
	place = copy.place ;
	isPresent = copy.isPresent ;
	successor = copy.successor ;
	pred = copy.pred ;
	route = copy.route ;
	params = copy.params ;
	moves = copy.moves ;
}

Noeud& Noeud::operator=(Noeud const& copy)
{
	// Copy constructor
	isDepot = copy.isDepot ;
	customer_or_depot_idx = copy.customer_or_depot_idx ;
	place = copy.place ;
	isPresent = copy.isPresent ;
	successor = copy.successor ;
	pred = copy.pred ;
	route = copy.route ;
	params = copy.params ;
	moves = copy.moves ;
	return *this;
}

Noeud::~Noeud(void){}

void Noeud::setRemaining()
{
	// seq1 has exactly the same meaning than seqi_j[0], but its more convenient to use and read
	seq1 = seqi_j[0] ;
	seq12 = seqi_j[1] ;
	seq123 = seqi_j[2] ;
	seq21 = seqj_i[1] ;
	seq321 = seqj_i[2] ;
}
