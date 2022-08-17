/**
 * @brief track the NOV of each node globally
 * @author Alexander Sanfilippo
 * @date 2022-08-02
 * @TODO test this function and update Genome.h/testing_grounds.cc to use it
 *     fix header hex code
 *     */
	

#ifndef NOV_H_87AD8052
#define NOV_H_87AD8052

#include <stdio.h>
#include <iostream>

#include <random>
#include <cstdlib>
//#include "gene.h"
//#include "Genome.h"
//#include "Species.h"
#include <unistd.h>
#include <vector>


std::default_random_engine gen(1996); //random values

using namespace std;
class NOV{
	private:
	vector <double> nov_vec = {0}; //stores each nodes nodal order values

	

	public:
	//CONSTRUCTORS
	/**
 	* @brief constructor that sets input and output nodes' NOVs*/
	NOV(int initial_inputs, int initial_outputs){
		for(int i =0; i < initial_inputs; i++){
			nov_vec.push_back(-1);
		}
		for(int i = 0; i < initial_outputs; i++){
			nov_vec.push_back(2);
		}
	}	
	
	//MEMBER FUNCTIONS
	/**
 	* @brief prints the NOV to console*/
	void summary(){
		cout << "The Global NOV:\n";
		for(int i = 0; i < int(nov_vec.size()); i++){
			cout << "node #" << i <<" has NOV: "<< nov_vec[i] << endl;
		}
	}
	/**
 	* @returns a nodes NOV given that node's number
 	* */
	double getNOV(int node_num){
		if(node_num <= int(nov_vec.size())-1 ){
			return nov_vec[node_num];
		}
		else{ //this node does not have an NOV yet
			return -1;
		}
	}
	int size(){
		return nov_vec.size();
	}
	/**
 	* @brief adds a new node and random NOV to the vector
 	* */
	void addNode(){
		//generate a new random nov
		auto urd = std::uniform_real_distribution<>(0,1); 
		//double rand_nov = urd(gen);
		//add the nov to the list
		nov_vec.push_back(urd(gen));	
	}
	
		
	/**
 	* @brief adds a new node with set NOV 
 	* */
	void addNode(double nov){
		//generate a new random nov
		//auto urd = std::uniform_real_distribution<>(0,1); 
		//double rand_nov = urd(gen);
		//add the nov to the list
		nov_vec.push_back(nov);	
	}

	
	/**
 	* @brief adds a new node with set NOV 
 	* */
	void addNode(int direction, double middle_nov){
		double nov;
		if(direction == -1){
			auto urd = std::uniform_real_distribution<>(0,middle_nov);
			nov = urd(gen);
		} 
		else{
			auto urd = std::uniform_real_distribution<>(middle_nov,1);
			nov = urd(gen);
		}
		nov_vec.push_back(nov);	
	}



	/**
 	* @brief check whether in/out node pair have correct order for cnx adding mutation
 	* */
	bool compareNodes(int in, int out){
		if(nov_vec[in] < nov_vec[out]){
			//correct order of nodes
			return true;		
		}
		else{
			return false;
		}
	}

	/**
 	* @brief averages the values of two nodes for node-adding mutation*/
	double avgNodes(int nodeA, int nodeB){
		double avg = (getNOV(nodeA) + getNOV(nodeB) ) / 2;
		return avg;

	}


};



#endif

