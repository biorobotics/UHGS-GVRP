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


#ifndef CLIENT_H
#define CLIENT_H

#include <string>
#include <vector>
#include <iostream>
#include "math.h"
using namespace std ;

// different client types
enum ClientType {AR_DEPOT, AR_CLIENT_NODE};

class Client
{
	public:

	// customer number
    int custNum ;

	// demand of a customer
	double demand ;

	// how many visits per week
	int freq ;

	// type of service
	ClientType ar_nodeType ;

	int ar_nbNodes ;
	vector <int> ar_Nodes ;

	// For a customer i, list of close customers j by proximity
	vector <int> closeCustomers ;

	// For a customer j, list of customers i to which j is considered as "close"
	vector <int> inverseCloseCustomers ;

	Client();

	~Client(void);
};

#endif
