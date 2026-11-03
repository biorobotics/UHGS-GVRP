//
// Created by chkwon on 3/22/22.
//

#ifndef INSTANCECVRPLIB_H
#define INSTANCECVRPLIB_H
#include<string>
#include<vector>
// #include <./../../spStruct.h>

class InstanceCVRPLIB
{
public:
	std::vector<double> x_coords;							//Cluster centroid coordinates
	std::vector<double> y_coords;							//Cluster centroid coordinates
	std::vector< std::vector<double> > dist_mtx;			//Minimum cluster distance matrix

	std::vector<double> service_time;
	std::vector<double> demands;
	double durationLimit = 1.e30;							// Route duration limit
	double vehicleCapacity = 1.e30;							// Capacity limit
	bool isDurationConstraint = false;						// Indicates if the problem includes duration constraints
	int nbClients ;											// Number of clients (excluding the depot)
	std::vector<std::vector<double>> arcs;					//Arcs of the full graph
	int depot;												//The depot is set to be the last node, so this is also the dimension
	std::vector<std::vector<std::vector<std::vector<double>>>> arcscluster; //Arcs between the clusters
	InstanceCVRPLIB(std::string pathToInstance, bool isRoundingInteger,bool cluvrp, bool isGVRP);
	int real_nbVeh;
	int n_nodes;
	bool cluvrp;											//Indicates if the file is a CluVRP instance with all the Hamiltonian paths or a standad GVRP
	std::vector<std::vector<double>> nodesCluLevel;			//Nodes obtained when transforming a CluVRP into a GVRP
	std::vector<std::vector<int>> cluster;					//List of nodes for each cluster
	std::vector < std::vector< double > > dist_mtx_nodes;	//Distance matrix node level
	std::vector< std::vector<double> > dist_mtxN;			//Distance matrix between the real nodes on the GVRP graph. It corresponds to dist_mtx_nodes if the 
															//instance is originally a GVRP
	std::vector<double> x_coordsN;							//Coordinates node level
	std::vector<double> y_coordsN;							//Coordinates node level
};


#endif //INSTANCECVRPLIB_H
