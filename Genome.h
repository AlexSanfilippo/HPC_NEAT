/**
 * @author Alex Sanfilippo
 * @date 2022-06-13
 * @breif define the Genome class
 *  This class holds a vector of genes which defines an
 *  individual, its network topology, and its neurons' weights
 * */

#ifndef GENOME_H_2376F5AB
#define GENOME_H_2376F5AB



#include <vector>
#include <iostream>
#include <fstream> //for writing to files
#include "gene.h"
#include "NOV.h"

//std::default_random_engine gen(1); //define random number engine "gen" //4321

class Genome{
	public:
		/*member variables*/
		std::vector <Gene> gene_vec; //holds the genes in a std vector
		double fitness = 1; //how fit this genome is on the test
		double adjusted = 1; //adjusted fitness score	


		/*Constructors, Destructors, etc...*/
		
		/**
 		* @brief creates efault initial network Genome
 		*
 		* */
		Genome(){
			//std::cout << "New Default Genome constructed\n";
			/*these four genes form the minimal network*/
			Gene A(1, 5, 1, true);
			A.setInnov(1);
			A.setInnovCount(1);
			Gene B(2, 5, 1, true);
			Gene C(3, 5, 1, true);
			Gene D(4, 5, 1, true);
			
			std::vector<Gene> temp_gene_vec = {A, B, C, D};	
			gene_vec = temp_gene_vec; //deep copy since no pointers in Gene objs 
			temp_gene_vec.clear();  //remove elements of vector
			temp_gene_vec.shrink_to_fit(); //deallocate memory of vector
			
			//std::cout << "InnovCount A " << A.getInnovCount() << "\n";
			//std::cout << "InnovCount B " << B.getInnovCount() << "\n";		
		}
		

	
		/**
 		* @brief constructs an EMPTY genome with NO genes
 		* */
		Genome(bool empty){
			//std::cout << "make an empty genome...\n";
		}

		/**
 		* @brief destructor
 		* */
		~Genome(){
			//std::cout << "Genome desctructor called...\n";
			//gene_vec.~vector(); //delete the gene_vec from memory --MEM ERROR
			
			//try clear and shrink_to_fit instead.  
			gene_vec.clear();
			gene_vec.shrink_to_fit();
		}
	
		/*Copy Constructor*/
		Genome(const Genome& old_g){
			//std::cout << "Genome copy constructor called\n";
			//this->gene_vec = old_g.gene_vec; //maybe causing segfaults
			for(int i = 0; i < old_g.size(); i++){
				gene_vec.push_back(old_g.gene_vec[i]);
			}
			//int size = old_g.size();
			//std::copy(&old_g.gene_vec, &old_g.gene_vec + size, &gene_vec);
			this->fitness = old_g.fitness;
			this->adjusted = old_g.adjusted;
		}
		/*copy Constructor old
		Genome(const Genome& old_g){
			std::cout << "Genome copy asgnmt Op called\n";
			gene_vec.clear(); //del the old gene_vec
			gene_vec = old_g.gene_vec;
			this->fitness = old_g.fitness;
                        this->adjusted = old_g.adjusted;
			return *this;
		}*/
		/*Move assignment operator*/
		Genome& operator=(Genome&& old_g){
			//std::cout << "Genome MV assmt op called\n";
			if (this != &old_g){
				//std::copy(&old_g.gene_vec, &old_g.gene_vec + old_g.size(), &gene_vec);
				//std::cout << "[Genome mv assmt op] after move\n";
				this->fitness = old_g.fitness;
				this->adjusted = old_g.adjusted;
						
				
				//neccessary to clear the left side's gene vec
				gene_vec.clear();
				gene_vec.shrink_to_fit();

				//iteratively deep copies old genes to new
				for(int i = 0; i < old_g.size(); i++){
					gene_vec.push_back(old_g.gene_vec[i]);
				}

				//std::cout << "[Genome mv assmt op] just before return...\n";
				//return *this;
			}
			//std::cout << "genome mv assmt op called on 2 identical genomes\n";
			//return *this;
			return *this;
		}
		

		//NEW cpy assmnt operator
		Genome& operator=(const Genome& other){
			//std::cout << "[Genome] copy assmt operator called...\n";
			if(this != &other){
				fitness = other.fitness;
				adjusted = other.adjusted;
				//gene_vec.swap(other.gene_vec); //is this okay? NO
				//std::copy(&other.gene_vec, &other.gene_vec + other.size(), &gene_vec); //bad alloc
				
		
				//neccessary to clear the left side's gene vec
				//std::cout << "attempting to clear and shrink this gene_vec\n";
				gene_vec.clear();
				gene_vec.shrink_to_fit();

				//std::cout << "attempting to push back genes\n";
				//iteratively deep copies old genes to new
				//std::cout << "call other.size()\n";
				//auto s = other.size();
				//std::cout << "other has size " << s << std::endl;
				for(int i = 0; i < other.size(); i++){
					gene_vec.push_back(other.gene_vec[i]);
				}
			}
			//std::cout << "returning *this\n";
			return *this;
		}
		/*tool functions*/
		int size() const { return gene_vec.size();}
		
		void get_gene_at(int index){
			Gene temp =  gene_vec[index];
			temp.summary();
		}
		
		void setFitness(double new_fitness){
			fitness = new_fitness;
		}
		double getFitness() const {
			return fitness;
		}
		double getAdjusted(){ return adjusted; };
		/**
 		* @brief prints out the genes and their weights/connections
 		*
 		* */
		void summary(){
			for(unsigned int i = 0; i < gene_vec.size(); i++){
				get_gene_at(i);
				std::cout << "\n";
			}
			std::cout << "adjusted = " << adjusted << std::endl;
		}
		
		/**
 		* @brief calculates and setsdjusted fitness score 
 		* */
		void calcAdjusted(int spec_pop){
			adjusted = fitness/spec_pop;	
		}
		
		/*dummyFitness function to use before we have neural nets set up*/
		void dummyFitness(){
			fitness = gene_vec.size();
		}
		
		/*
 		* @brief wrtes the node_inputs and node_weights vectors to a file for keras to read
 		*
 		* **/
		void  writeToFile(ofstream *fp, int genome_index){
			
			//*fp << "START" << "\n";
			

			//write the genome's nodal info to the .csv file	
			for(int i = 1; i < int( node_inputs.size() ); i++){ //for each node
				if(node_inputs[i][0] > 0){//if node has at least 1 input
					*fp << i << ","; //write node number
					//print the input-nodes' numbers
					for(int j = 0; j < int(node_weights[i].size()+1); j++){ //for each nodes inputs
						*fp << node_inputs[i][j] << ","; //write input-node number 
					}
					//print the weights	
					for(int j = 0; j < int(node_weights[i].size()); j++){ //for each nodes inputs
						*fp << node_weights[i][j] << ","; //write weight
					}
					*fp << "\n";
				}
			}	
			*fp << "END," << "\n";

	
		}
		/**
 		* @brief makes a vector of all the Genome's node's and their weights
 		* */
		void calcNodeInfo(NOV *nov){
			const int highest_node = int(nov->size());
			for(int i = 0; i < highest_node; i++){
				std::vector <double> row = {0};
				node_inputs.push_back(row);
				std::vector <double> row_w;
                                node_weights.push_back(row_w);
			}
			for(int i = 0; i < int(gene_vec.size()); i++){ //for each gene
				int input = gene_vec[i].getIn();
				int output = gene_vec[i].getOut();
				double w = gene_vec[i].getWeight();
				bool enabled = gene_vec[i].getEnabled();
				
				if(enabled){ //only enabled genes are built on the NN
					node_inputs[output][0] += 1; //increment the input-node count
					node_inputs[output].push_back(input);  //add this input-node to the list
					node_weights[output].push_back(w); //add weight to the weight list
				}
			}
			//std::cout << "finished populating the vectors\n";
			
			/*
			//TP print the ni and nw vecs to check this is correct 
			for(int i = 1; i < highest_node; i++){ //for each node
				std::cout << "node #" << i << ": "; 
				//print the input-nodes' numbers
				for(int j = 0; j < int(node_weights[i].size()+1); j++){ //for each nodes inputs
					std::cout << node_inputs[i][j] << ", "; 
				}
				//print the weights	
				for(int j = 0; j < int(node_weights[i].size()); j++){ //for each nodes inputs
					std::cout << node_weights[i][j] << ", "; 
				}
				std::cout << "\n";
			}	
			*/

		}
		/**
 		* @brief writes the Genome to a file for python to use
 		* */
		/*
		void writeToFile(fp){
			
		}
		*/

		/***
 		* @brief call to mutate the genes
 		* @param rng_vec a pointer to a vector of random doubles
 		* */
		void mutate(NOV *nov){
			//std::cout << "inside mutate()\n";
			/*WIP: IN-Function RNG generation*/
			//auto uid = std::uniform_int_distribution<>(0,10); //value between 0 and 10, inclusive?
			//uid(gen);
			//
			//distribution for normal, doubles, 0 to 1 (probabilities)
			auto urd = std::uniform_real_distribution<>(0,1); //random probability
			auto uid = std::uniform_int_distribution<>(0,gene_vec.size() - 1); //random gene
			/*NOTE: may need min to be the highest node index of the input nodes, not sure if
 			* this will cause problems otherwise.*/
			//auto uid_node = std::uniform_int_distribution<>(1, nov->size()); //random node(neuron)
			//auto uid_node = std::uniform_int_distribution<>(1, nov->size()); //random INPUT node(neuron)
			auto urd_weight = std::uniform_real_distribution<>(-MUT_MAGN, MUT_MAGN); 
			auto urd_p = std::uniform_real_distribution<>(-2,2); //for weight perturbation
			if(urd(gen) < P_weight){ //if we mutate this GENOME's weights
				
				/*CHANGE WEIGHTS OF NODES*/
				for(unsigned long int i = 0; i < gene_vec.size(); i++){ //for each gene
					if(urd(gen) < P_perturb){ //perturb old weight
						double old_weight = gene_vec[i].getWeight();
                                                gene_vec[i].setWeight(old_weight + old_weight*urd_p(gen));	
						
					/*/std::cout << "PERTURBED gene #" << gene_vec[i].getInnov() << "'s weight to " \
                                                << gene_vec[i].getWeight() << std::endl;  */
					}
					else{ //new weight
						gene_vec[i].setWeight(urd_weight(gen));
						/*
						std::cout << "SET gene #" << gene_vec[i].getInnov() <<\ "'s weight to "\						<< gene_vec[i].getWeight() << std::endl; */
					}
				}
			}

			/*ADD NEW NODE*/
			
			if(urd(gen) < P_newnode){ //if we create a new node
				int cur_index = uid(gen);
			//	int cur_innov = gene_vec[cur_index].getInnov();
				//std::cout << "creating new node @ gene #" << cur_innov << std::endl;
				//split the cur_index gene into 2 new genes, removing the original
						int in = gene_vec[cur_index].getIn();
						int out = gene_vec[cur_index].getOut();
						double weight = gene_vec[cur_index].getWeight();
						//node_count += 1; //increment count of nodes defined by Genome
						node_count = nov->size();
						Gene gene_a(in, node_count , 1, true); 
						//std::cout << "ADDNODE mutation added cnx A tween " << in <<"-->"<<node_count<<"\n";
						Gene gene_b(node_count, out, weight, true );
						//std::cout << "ADDNODE mutation added cnx B tween " << node_count <<"-->"<<out<<"\n";

						/*std::cout << "created 2 new genes with in #'s" << gene_a.getInnov() <<\
							gene_b.getInnov() << std::endl;
						*/
						//disable the old node
						gene_vec[cur_index].setEnabled(false);

						//delete the old node from the Genome	--WRONG
						//gene_vec.erase(gene_vec.begin() + cur_index);
						//gene_vec.shrink_to_fit();

						//add new vectors to the genome	
						gene_vec.push_back(gene_a);
						gene_vec.push_back(gene_b);	


						/*Add New Node to NOV vector*/
						//std::cout << "in = " << in << ", out = " << out << endl;
						if(nov->getNOV(in) != -1 && nov->getNOV(out) != 2){ 
						//if neither in or out node is a final input or output
							nov->addNode(nov->avgNodes(in,out)); //add nov to the NOV object's vector
							//std::cout << "adding node tween 2 hidden nodes\n"; 
						}
						else if(nov->getNOV(in) == -1 && nov->getNOV(out) == 2){
							//std::cout << "adding node tween 2 final in/out layers\n";
							nov->addNode(); //add a random nov 
						}
						else if(nov->getNOV(in) == -1 && nov->getNOV(out) != 2){ //input to hidden
							//std::cout << "add value tween input and hidden\n";
							nov->addNode(-1, nov->getNOV(out));
						}
						else if(nov->getNOV(in) != -1 && nov->getNOV(out) == 2){ //hidden to output
							//std::cout << "add node tween hidden and output\n";
							nov->addNode(2, nov->getNOV(in));
						}
			}
			
			/*ADD NEW LINK (CNXN)*/
			
			if(urd(gen) < P_newlink){
				
				/* OLD RANDOM INPUT/OUPUT node picker- results in floating/ghost nodes
 					* ie: nodes with no inputs, only outputs.  remove to up
 					* efficiency	 		
				int rand_in = 0;
				int rand_out = 0;
				//if in/out pair are the same value twice or BOTH INPUTS, redraw
				while(rand_in == rand_out || (rand_in <= 4 && rand_out <= 4) ){
					std::cout << rand_in<<"=="<<rand_out << ", re-rolling values\n";
					rand_in = uid_node(gen);  
					rand_out = uid_node(gen);
				}
				*/
				
				//NEW RANDOM INPUT/OUPUT node picker
				
				int cur_indexA = 0; //draw a random gene
				int cur_indexB = 0;
				while(cur_indexA == cur_indexB){
					cur_indexA = uid(gen);
					cur_indexB = uid(gen);
				}
				int rand_in = 0;
				int rand_out = 0; 
				//randomly decide whether to use gene A's input
				//or output as the new cnxn's input. same
				//logic for output from gene "B"
				while(rand_in == rand_out || (rand_in <= 4 && rand_out <= 4) ){
					if(urd(gen) < 0.5){
						rand_in = gene_vec[cur_indexA].getIn();
					}	
					else{
						rand_in = gene_vec[cur_indexA].getOut();
					}
					if(urd(gen) < 0.5){
						rand_out = gene_vec[cur_indexB].getIn();
					}
					else{
						rand_out = gene_vec[cur_indexB].getOut();
					}
				}
				/*Check that NOV values in right order, otherwise reverse in and out*/
				if(nov->getNOV(rand_in) > nov->getNOV(rand_out)){
					//std::cout << "re-ordering rand_in and rand_out\n";
					int temp = rand_in;
					rand_in = rand_out;
					rand_out = temp;
					//std::cout << rand_in<<" =in, out = "<<rand_out << "\n";
	
				}

				/*Check if the new in-out pair already exists in a gene in this Genome*/
				bool do_create = true;
				//loop thru the gene vec
				for(int i = 0; i < int(gene_vec.size()); i++){
					//if the input and output are the same
					if(rand_out == gene_vec[i].getOut() && rand_in == gene_vec[i].getIn()){
						//cancel this mutation
						do_create = false;
						break;
					}
				}
				
				

				if(do_create){	
					//create new Gene for this link
					//Gene gene_new(rand_in, rand_out, urd_weight(gen), true);
					Gene gene_new(rand_in, rand_out, 0, true);
					//std::cout << "NEWLINK mutation created link " << rand_in <<"-->"<<rand_out<<std::endl;
					gene_vec.push_back(gene_new);	
					/*std::cout << "added new link with in #'s" << gene_new.getInnov() <<\
								std::endl;
					*/
					//std::cout << "established new link between" << rand_in << " and " << rand_out << std::endl;
				}		
				



				
			}
		}

			




	private:
		
		//Contain Genome's nodal data for keras to use.  write this to file
		std::vector <std::vector <double>> node_inputs; //all nodes and thier input nodes
		std::vector <std::vector<double>> node_weights; //all nodes and all thier weights

		//std::vector <Gene> gene_vec; //holds the genes in a std vector 
	//	double fitness; //how fit this genome is on the test
	//	double adjusted; //adjusted fitness score	
		const int MUT_MAGN = 100;  /*max ammount weight can change*/
		int node_count = 5; //be sure to set according to initial architecture 
		
		/*MUTATION PROBABILITY THRESHOLDS*/
		double P_weight = 0.8; //chance this genome will mutate its weight(s)
		double P_perturb = 0.9; //chance to change a genomes weight by a multiple 
		double P_newweight = 0.1; //chance to replace weight with new random weight
		double P_newnode = 0.05; //add new gene for new node, splitting one into two //should be 0.05
		double P_newlink = 0.30; //add new connection (link) between 2 existing nodes (new gene) 
			
		/*Other Hypers*/
};


#endif


