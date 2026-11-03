//
// Created by chkwon on 3/22/22.
//

#include <fstream>
#include <cmath>
#include "InstanceCVRPLIB.h"


#include <iostream>
#include <vector>
#include <queue>
#include <climits>
#include <algorithm>

InstanceCVRPLIB::InstanceCVRPLIB(std::string pathToInstance, bool isRoundingInteger = true, bool cluvrp = false, bool isGVRP=false)
{
		std::string content, content2, content3;
		double serviceTimeData = 0.;
		//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
		//isDurationConstraint = true;
		//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
		// Read INPUT dataset
		std::ifstream inputFile(pathToInstance);
		if (inputFile.is_open())
		{
			getline(inputFile, content);
			getline(inputFile, content);
			//getline(inputFile, content);
			for (inputFile >> content; content != "NODE_COORD_SECTION"; inputFile >> content)
			{
				if (content == "GVRP_SETS") { inputFile >> content2 >> nbClients; } // Need to substract the depot from the number of nodes
				else if (content == "EDGE_WEIGHT_TYPE")	inputFile >> content2 >> content3;
				else if (content == "CAPACITY")	inputFile >> content2 >> vehicleCapacity;
				else if (content == "DISTANCE") { inputFile >> content2 >> durationLimit; isDurationConstraint = true; }
				else if (content == "SERVICE_TIME")	inputFile >> content2 >> serviceTimeData;
				else if (content == "VEHICLES") inputFile >> content2 >> real_nbVeh;
				else if (content == "DIMENSION") inputFile >> content2 >> n_nodes;

				else throw std::runtime_error("Unexpected data in input file: " + content);
			}
			if (nbClients <= 0) throw std::runtime_error("Number of nodes is undefined");
			if (vehicleCapacity == 1.e30) throw std::runtime_error("Vehicle capacity is undefined");

			x_coords = std::vector<double>(n_nodes);
			y_coords = std::vector<double>(n_nodes);
			demands = std::vector<double>(nbClients + 1);
			service_time = std::vector<double>(nbClients + 1);
			cluster = std::vector<std::vector<int>>(nbClients+1);
			depot = n_nodes + 1;
			double eps = 0.00001;
			// Reading node coordinates
			// depot must be the first element
			// 		- i = 0 in the for-loop below, or
			// 		- node_number = 1 in the .vrp file
			// customers are
			// 		- i = 1, 2, ..., nbClients in the for-loop below, or
			// 		- node_number = 2, 3, ..., nb_Clients in the .vrp file
			int node_number;
			for (int i = 0; i < n_nodes; i++)
			{
				inputFile >> node_number >> x_coords[i] >> y_coords[i];
				if (node_number != i + 1) throw std::runtime_error("The node numbering is not in order.");
			}
			// Reading cluster membership information
			inputFile >> content;
			if (content != "GVRP_SET_SECTION") throw std::runtime_error("Unexpected data in input file: " + content);
			int i = 1;
			int cluster_number;
			inputFile >> cluster_number;
			while( i< nbClients+1)
			{
				int value;
				while (inputFile >> value)
				{
					if (value == -1) {
						i++;
						if (i == nbClients + 1) {
							break;
						}

						inputFile >> cluster_number;
						if (cluster_number != i ) throw std::runtime_error("The cluster numbering is not in order.");

						break;
					}
					cluster[i].push_back(value-1);
				}
			}

			std::vector <int> nodes_cluster= std::vector <int>(n_nodes);
			nodes_cluster[0] = 0;
			for (int i = 0; i < nbClients+1; i++)
			{
				for (int j = 0; j < cluster[i].size(); j++)
				{
					nodes_cluster[cluster[i][j]] = i;
				}
			}



			// Reading demand information
			inputFile >> content;
			if (content != "DEMAND_SECTION") throw std::runtime_error("Unexpected data in input file: " + content);
			demands[0] = 0;
			for (int i = 1; i <= nbClients; i++)
			{
				inputFile >> content >> demands[i];
				service_time[i] = (i == 0) ? 0. : serviceTimeData;
			}
			arcscluster = std::vector<std::vector<std::vector<std::vector<double>>>>(nbClients + 1, std::vector<std::vector<std::vector<double>>>(nbClients + 1));
			std::vector<double>x_coords_clu(nbClients + 1, 0);
			std::vector<double>y_coords_clu(nbClients + 1, 0);
			std::vector<int>n_clu(nbClients + 1, 0);

			if (cluvrp) {
				std::vector < std::vector< double > > nodesCluLevel1(n_nodes, std::vector<double>(n_nodes));
				//Reading intra cluster distances 
				nodesCluLevel.resize(n_nodes * n_nodes, std::vector<double>(4));
				inputFile >> content;
				if (content != "INTRA_CLUSTER_DISTANCE" && content != "EOF") throw std::runtime_error("Unexpected data in input file: " + content);
				int k = 1;
				for (int i = 0; i < n_nodes; i++)
				{
					for (int j = 0; j < n_nodes; j++)
					{
						double value;
						inputFile >> value;
						if (value > 10000) {
							continue;
						}
						nodesCluLevel[k]={ (double)i,(double)j, value,(double)nodes_cluster[i] };
						k++;
						nodesCluLevel1[i][j] = value;
					}
				}

				for (int i = 1; i < nbClients+1; i++)
				{
					if (cluster[i].size() < 2) {
						nodesCluLevel[k] = { (double)cluster[i][0],(double)cluster[i][0], 0,(double)i };
						k++;
					}
				}


				nodesCluLevel.erase(nodesCluLevel.begin() + k, nodesCluLevel.end());

				int n_nodesC = nodesCluLevel.size()+1;
				depot = n_nodesC + 1;

				// Calculating the distance matrix

				dist_mtx = std::vector < std::vector< double > >(nbClients + 1, std::vector <double>(nbClients + 1));
				dist_mtx_nodes = std::vector < std::vector< double > >(n_nodes, std::vector <double>(n_nodes));

				for (int i = 0; i < n_nodes; i++)
				{
					for (int j = 0; j < n_nodes; j++)
					{
						dist_mtx_nodes[i][j] = std::sqrt(
							(x_coords[i] - x_coords[j]) * (x_coords[i] - x_coords[j])
							+ (y_coords[i] - y_coords[j]) * (y_coords[i] - y_coords[j])
						);

						if (isRoundingInteger) dist_mtx_nodes[i][j] = round(dist_mtx_nodes[i][j]);
					}
				}

				//Mark the arcs that should be used
				std::vector < std::vector< int > > marked(n_nodes, std::vector<int>(n_nodes));
				for (int i = 1; i < nbClients+1; i++)
				{

					for (int j = 1; j < nbClients+1; j++)
					{
						if (i == j) continue;

						if (cluster[i].size() < 3)
						{
							for (int p = 0; p < cluster[i].size(); p++) {
								for (int q = 0; q < cluster[j].size(); q++) {

									marked[cluster[i][p]][cluster[j][q]]=1;
								}
							}
							continue;
						}

						if (cluster[j].size() < 3)
						{
							for (int p = 0; p < cluster[i].size(); p++) {
								for (int q = 0; q < cluster[j].size(); q++) {

									marked[cluster[i][p]][cluster[j][q]] = 1;
								}
							}
							continue;
						}

						for (int p = 0; p < cluster[i].size(); p++)
						{
							int idx_p = cluster[i][p];
							for (int q = 0; q < cluster[j].size(); q++)
							{
								int idx_q = cluster[j][q];
								double alpha = 1.e30;
								std::vector < std::vector< double > > alpha_vect(cluster[j].size() * cluster[i].size(), std::vector<double>(3));
								int k_al = 0;
								for (int ii = 0; ii < cluster[i].size(); ii++)
								{
									if (ii == p) continue;
									int idx_i = cluster[i][ii];
									for (int jj = 0; jj < cluster[j].size(); jj++)
									{
										if (jj == q) continue;
										int idx_j = cluster[j][jj];

										//double c_ij = nodesCluLevel1[idx_i][idx_j];
										alpha_vect[k_al] = { (double)idx_i,(double)idx_j,dist_mtx_nodes[idx_i][idx_j] + nodesCluLevel1[idx_p][idx_i] + nodesCluLevel1[idx_j][idx_q] };
										k_al++;
									}
								}
								alpha_vect.erase(alpha_vect.begin() + k_al, alpha_vect.end());
								for (int iii = 0; iii < alpha_vect.size(); iii++)
								{
									if (alpha_vect[iii][2] < alpha) {
										alpha = alpha_vect[iii][2];
									}
								}
								for (int iii = 0; iii < alpha_vect.size(); iii++)
								{
									if (alpha_vect[iii][2] - alpha<=eps) {
										marked[alpha_vect[iii][0]][alpha_vect[iii][1]]=1;
									}
								}
							}
						}
					}
				}

				//arcs between clusters
				std::vector < std::vector< double > > arcsCluLevel(n_nodesC* n_nodesC, std::vector< double >(5));
				k = 0;
				for (int i = 0; i < nodesCluLevel.size(); i++)
				{
					for (int j = i + 1; j < nodesCluLevel.size(); j++)
					{
						if (nodesCluLevel[i][3] == nodesCluLevel[j][3]) continue;
						if (marked[nodesCluLevel[i][1]][nodesCluLevel[j][0]] > 0) {
							arcsCluLevel[k] = { (double)i , (double)j ,nodesCluLevel[i][2] + dist_mtx_nodes[nodesCluLevel[i][1]][nodesCluLevel[j][0]],nodesCluLevel[i][3] ,nodesCluLevel[j][3] };
							k++;
						}
						if (marked[nodesCluLevel[j][1]][nodesCluLevel[i][0]] > 0) {
							arcsCluLevel[k] = { (double)j , (double)i ,nodesCluLevel[j][2] + dist_mtx_nodes[nodesCluLevel[j][1]][nodesCluLevel[i][0]],nodesCluLevel[j][3] ,nodesCluLevel[i][3] };
							k++;
						}
					}
				}
				//arcs from depot to cluster
				for (int i = 0; i < nodesCluLevel.size(); i++)
				{
					arcsCluLevel[k]={ (double)i , 0,nodesCluLevel[i][2] + dist_mtx_nodes[nodesCluLevel[i][1]][0], nodesCluLevel[i][3] , 0 };
					k++;
					arcsCluLevel[k]={ 0, (double)i , dist_mtx_nodes[0][nodesCluLevel[i][0]], 0 , nodesCluLevel[i][3] };
					k++;
				}
				arcsCluLevel.erase(arcsCluLevel.begin() + k, arcsCluLevel.end());



				// Computing the arcscluster tensor
				for (int i = 0; i < arcsCluLevel.size(); i++)
				{
					if (arcsCluLevel[i][0] == 2685){
						int ii=4;
					}
					arcscluster[arcsCluLevel[i][3]][arcsCluLevel[i][4]].push_back({ static_cast<double>(arcsCluLevel[i][0]) ,static_cast<double>(arcsCluLevel[i][1]) ,arcsCluLevel[i][2] });
				}
				
				// Computing the distance matrix cluster level
				dist_mtxN.resize(n_nodesC, std::vector<double>(n_nodesC,1.e10));
				for (int i = 0; i < arcsCluLevel.size(); i++)
				{
					dist_mtxN[static_cast<double>(arcsCluLevel[i][0])][static_cast<double>(arcsCluLevel[i][1])] = arcsCluLevel[i][2];
				}


				for (int i = 0; i < nbClients + 1; i++)
				{
					for (int j = 0; j < nbClients + 1; j++)
					{
						if (i == j) {
							dist_mtx[i][j] = 0;
							continue;
						}
						double incumbent = 1.e30;
						for (int k = 0; k < arcscluster[i][j].size(); k++) {
							if (arcscluster[i][j][k][2] < incumbent) {
								incumbent = arcscluster[i][j][k][2];
							}
						}
						dist_mtx[i][j] = incumbent;
					}
				}


				//Compute the coodinates of the centroids for clusters
				

				for (int i = 0; i < n_nodes; i++)
				{
					x_coords_clu[nodes_cluster[i]] += x_coords[i];
					y_coords_clu[nodes_cluster[i]] += y_coords[i];
					n_clu[nodes_cluster[i]] = n_clu[nodes_cluster[i]] + 1;
				}
				for (int i = 0; i < nbClients; i++)
				{
					x_coords_clu[i] = x_coords_clu[i] / n_clu[i];
					y_coords_clu[i] = y_coords_clu[i] / n_clu[i];

				}

				x_coordsN = x_coords;
				y_coordsN = y_coords;

				x_coords = x_coords_clu;
				y_coords = y_coords_clu;




				// Reading depot information (in all current instances the depot is represented as node 1, the program will return an error otherwise)
				inputFile >> content >> content2 >> content3 >> content3;
				//			if (content != "DEPOT_SECTION") throw std::runtime_error("Unexpected data in input file: " + content);
				//			if (content2 != "1") throw std::runtime_error("Expected depot index 1 instead of " + content2);
				//if (content != "EOF") throw std::runtime_error("Unexpected data in input file: " + content3);




			}






			//VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
			else {//If it is not a golden/CluVRP instance
				// Calculating 2D Euclidean Distance
				dist_mtx = std::vector < std::vector< double > >(nbClients + 1, std::vector <double>(nbClients + 1));
				std::vector < std::vector< double > > dist_mtx_nodes = std::vector < std::vector< double > >(n_nodes, std::vector <double>(n_nodes));

				for (int i = 0; i < n_nodes; i++)
				{
					for (int j = 0; j < n_nodes; j++)
					{
						dist_mtx_nodes[i][j] = std::sqrt(
							(x_coords[i] - x_coords[j]) * (x_coords[i] - x_coords[j])
							+ (y_coords[i] - y_coords[j]) * (y_coords[i] - y_coords[j])
						);

						if (isRoundingInteger) dist_mtx_nodes[i][j] = round(dist_mtx_nodes[i][j]);
					}
				}

				//Mark the arcs that should be used
				std::vector < std::vector< int > > marked(n_nodes, std::vector<int>(n_nodes));
				for (int i = 1; i < nbClients + 1; i++)
				{

					for (int j = 1; j < nbClients + 1; j++)
					{
						if (i == j) continue;

						if (cluster[i].size() < 3)
						{
							for (int p = 0; p < cluster[i].size(); p++) {
								for (int q = 0; q < cluster[j].size(); q++) {

									marked[cluster[i][p]][cluster[j][q]] = 1;
								}
							}
							continue;
						}

						if (cluster[j].size() < 3)
						{
							for (int p = 0; p < cluster[i].size(); p++) {
								for (int q = 0; q < cluster[j].size(); q++) {

									marked[cluster[i][p]][cluster[j][q]] = 1;
								}
							}
							continue;
						}



						for (int p = 0; p < cluster[i].size(); p++)
						{
							int idx_p = cluster[i][p];
							for (int q = 0; q < cluster[j].size(); q++)
							{
								int idx_q = cluster[j][q];
								double alpha = 1.e30;
								std::vector < std::vector< double > > alpha_vect(cluster[j].size() * cluster[i].size(), std::vector<double>(3));
								int k_al = 0;
								for (int ii = 0; ii < cluster[i].size(); ii++)
								{
									if (ii == p) continue;
									int idx_i = cluster[i][ii];
									for (int jj = 0; jj < cluster[j].size(); jj++)
									{
										if (jj == q) continue;
										int idx_j = cluster[j][jj];

										alpha_vect[k_al] = { (double)idx_i,(double)idx_j,dist_mtx_nodes[idx_i][idx_j] + dist_mtx_nodes[idx_p][idx_i] + dist_mtx_nodes[idx_j][idx_q] };
										k_al++;
									}
								}
								alpha_vect.erase(alpha_vect.begin() + k_al, alpha_vect.end());
								for (int iii = 0; iii < alpha_vect.size(); iii++)
								{
									if (alpha_vect[iii][2] < alpha) {
										alpha = alpha_vect[iii][2];
									}
								}
								for (int iii = 0; iii < alpha_vect.size(); iii++)
								{
									if (alpha_vect[iii][2] - alpha <= eps) {
										marked[alpha_vect[iii][0]][alpha_vect[iii][1]] = 1;
									}
								}
							}
						}
					}
				}

				for (int i = 0; i < n_nodes; ++i) {
					marked[0][i] = 1;
					marked[i][0] = 1;
				}


				// Computing the arcscluster tensor
				//
				for (int i = 0; i < n_nodes; ++i) {
					for (int j = 0; j < n_nodes; j++)
					{
						if (nodes_cluster[i] == nodes_cluster[j]) continue;
						if (marked[i][j] > 0) {
							arcscluster[nodes_cluster[i]][nodes_cluster[j]].push_back({ static_cast<double>(i) ,static_cast<double>(j) ,dist_mtx_nodes[i][j] });
						}
					}

				}
				// Computing the distance matrix cluster level
				for (int i = 0; i < nbClients + 1; i++)
				{
					for (int j = 0; j < nbClients + 1; j++)
					{
						if (i == j) {
							dist_mtx[i][j] = 0;
							continue;
						}
						double incumbent = 1.e30;
						for (int k = 0; k < arcscluster[i][j].size(); k++) {
							if (arcscluster[i][j][k][2] < incumbent) {
								incumbent = arcscluster[i][j][k][2];
							}
							//	incumbent += arcscluster[i][j][k][2];
						}
						dist_mtx[i][j] = incumbent;
					}
				}
				dist_mtxN = dist_mtx_nodes;


				//Compute the coodinates of the centroids for clusters
				std::vector<double>x_coords_clu(nbClients + 1, 0);
				std::vector<double>y_coords_clu(nbClients + 1, 0);
				std::vector<int>n_clu(nbClients + 1, 0);


				for (int i = 0; i < n_nodes; i++)
				{
					x_coords_clu[nodes_cluster[i]] += x_coords[i];
					y_coords_clu[nodes_cluster[i]] += y_coords[i];
					n_clu[nodes_cluster[i]] = n_clu[nodes_cluster[i]] + 1;
				}
				for (int i = 0; i < nbClients; i++)
				{
					x_coords_clu[i] = x_coords_clu[i] / n_clu[i];
					y_coords_clu[i] = y_coords_clu[i] / n_clu[i];

				}

				x_coords = x_coords_clu;
				y_coords = y_coords_clu;

				// Reading depot information (in all current instances the depot is represented as node 1, the program will return an error otherwise)
				inputFile >> content >> content2 >> content3 >> content3;
				//			if (content != "DEPOT_SECTION") throw std::runtime_error("Unexpected data in input file: " + content);
				//			if (content2 != "1") throw std::runtime_error("Expected depot index 1 instead of " + content2);
				if (content != "EOF") throw std::runtime_error("Unexpected data in input file: " + content3);				// Calculating 2D Euclidean Distance
				dist_mtx = std::vector < std::vector< double > >(nbClients + 1, std::vector <double>(nbClients + 1));
				//dist_mtx_nodes = std::vector < std::vector< double > >(n_nodes, std::vector <double>(n_nodes));

				// Computing the arcscluster tensor
				for (int i = 0; i < n_nodes; ++i) {
					for (int j = 0; j < n_nodes; j++)
					{
						if (nodes_cluster[i] == nodes_cluster[j]) continue;
						arcscluster[nodes_cluster[i]][nodes_cluster[j]].push_back({ static_cast<double>(i) ,static_cast<double>(j) ,dist_mtx_nodes[i][j] });

					}

				}
				// Computing the distance matrix cluster level
				for (int i = 0; i < nbClients + 1; i++)
				{
					for (int j = 0; j < nbClients + 1; j++)
					{
						if (i == j) {
							dist_mtx[i][j] = 0;
							continue;
						}
						double incumbent = 1.e30;
						for (int k = 0; k < arcscluster[i][j].size(); k++) {
							if (arcscluster[i][j][k][2] < incumbent) {
								incumbent = arcscluster[i][j][k][2];
							}
						}
						dist_mtx[i][j] = incumbent;
					}
				}
				dist_mtxN = dist_mtx_nodes;



				// Reading depot information (in all current instances the depot is represented as node 1, the program will return an error otherwise)
				inputFile >> content >> content2 >> content3 >> content3;
				//			if (content != "DEPOT_SECTION") throw std::runtime_error("Unexpected data in input file: " + content);
				//			if (content2 != "1") throw std::runtime_error("Expected depot index 1 instead of " + content2);
				if (content != "EOF") throw std::runtime_error("Unexpected data in input file: " + content3);
			}
		}
		else
			throw std::runtime_error("Impossible to open instance file: " + pathToInstance);
	
	}
