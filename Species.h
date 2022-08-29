/** @author Alexander Sanfilippo
 *  * @date 11-06-2022
 *   * @TODO clean up the code, add additional comments
 */

#ifndef SPECIES_H_8f60af27
#define SPECIES_H_8f60af27

#include <vector>
#include <math.h>
#include "gene.h"
#include "Genome.h"


std::default_random_engine mygen(1996);

class Species {

        public:
		
		//Constructors
		
		Species(int innov_count){
			Gene g = Gene(0,0,0,true);
			g.setInnovCount(innov_count);
			name = "1";	
		}
		/*for making a new species given a genome that does not
 		* fit into any existing species.  Not sure if need a copy
 		* constructor for Genome for this to work right.  Issue 
 		* may arise with deep vs shallow copying.*/
		Species(Genome g, int innov_count, string parent_name, int parent_subspecies_count){	
			genome_vec.push_back(g);
			rep = g; 
			rep_index = 0;

			//update innov_count 
			Gene g2 = Gene(0,0,0,true);
                        g2.setInnovCount(innov_count);
	
			//give this species a name based on parents name
			name = parent_name + "_" + to_string(parent_subspecies_count+1);

		}

		std::vector <Genome> genome_vec; //vector of genomes in this species
		

		/**
 		* @brief returns the highest fitness achieved by this species*/
		double getHighestFitness(){
			return highest_fitness;
		} 
		
		/**
 		* @breif getter function for last_improved, 
 		* */
		int getLastImproved(){
			return last_improved;
		}
			
		/**
 		* @breif setter function for last_improved, 
 		* */
		void setLastImproved(int last){
			last_improved = last;
		}
		
		/**
 		* @breif setter function for last_improved, 
 		* */
		void setHighestFitness(double fit){
			highest_fitness = fit;
		}
		
		
		/**
 		 * @brief check if new species fitness record has been acheived
 		 * @param current_fitness the average fitness of this species this Gen.
 		 *
 		 */
		void checkImprovement(double current_fitness){
			if(current_fitness < highest_fitness){ //if no improvement
				last_improved += 1; //count number of gens. since last improvement
			}
			else{
				last_improved = 0; //reset generation counter
				highest_fitness = current_fitness; //set new best species fitness
			}
		}

					
		/**
 		* @brief returns this species name
 		* @return std::string name*/
		string get_name(){return name; }
		
		/**
 		* @brief sets this species name
 		* @return void*/
		void set_name(string n){
			name = n; 
		}
		/**
 		* @brief return number of subspecies that broke away from this species
 		* */
		int get_num_subspecies(){return num_subspecies; }
		
		/**
 		* @breif set number of subspecies
 		* */
		void set_num_subspecies(int n){
			num_subspecies = n;
		}
		/**
 		* @brief increment the number of subspecies
 		* */
		void add_subspecies(){
			num_subspecies += 1;
		}

		//DEATH CLOCK FXNS
		void set_death_clock(int n){
			death_clock = n;
		}	
		int get_death_clock(){
			return death_clock;
		}
		void add_death_clock(){
			death_clock += 1;
		}
		


		Genome rep; //copy of representative genome
		int rep_index; //index of representative genome


		/**
 		* @brief randomly chooses a new representative genome for speciation
 		*
 		* */
		void chooseRep(){
			auto uid = std::uniform_int_distribution<>(0,genome_vec.size() - 1); 
			int rand_index = uid(mygen);
			rep_index = rand_index;
			rep = genome_vec[rand_index]; 
		}
		/**
 		* @brief Wrapper for std::vector size function
 		* */
		int size(){
			return int(genome_vec.size());
		}
		/**
 		* @brief sets the avg fitness of the species
 		* */
		void setFitness(double f){
			fitness = f;
		}
		
		/**
 		* @brief gets the avg fitness of the species
 		* */
		double getFitness(){
			return fitness;
		}
		
		/**
 		* @brief calculates the avg fitness of the species
 		* */
		void calcFitness(){
			int size = this->size();
			if(size > 0){
				double avg = 0;
				for(int i = 0; i < size; i++){
					avg += genome_vec[i].getFitness();
				}
				avg = avg/size;
				fitness = avg;	
			}
			else{ //handles empty species
				fitness = 0;
			}
		}
		/**
 		* @brief adds genome g to the genome_vec of this species*/
		void addGenome(Genome g){
			genome_vec.push_back(g);

		}
		//Wrapper
		void push_back(Genome g){
			genome_vec.push_back(g);
		}
	
		/** 
 		* @brief determines sum of adjusted fitnesses for all genomes
 		* */
		double calcSumOfAdjusted(){
			//std::cout << "inside calcSumOfAdjusted() for genome with size " << genome_vec.size() <<std::endl;
			double sum = 0;
			for(unsigned long int i = 0; i < genome_vec.size(); i++){
				sum += genome_vec[i].adjusted;
			}

			return sum;
		}
		void setOffspring(int count){
			offspring = count;
		}
		int getOffspring() const { return offspring; };
		
		/**
 		* @brief sorts the Species by fitness
 		* uses quick sort algorithm - used wikipedia pseudo as guide
 		* */
		void sort(std::vector <Genome>& V, int low, int high){
			//std::cout << "inside sort() low=" <<low<<", high="<<high << std::endl;
			//std::vector <Genome> *Vptr = genome_vec; //create pointer V to the genome_vec of this species
			//std::vector <Genome> &V = *Vptr; //create a reference (for normal V[] access
			if( low >= 0 && high >= 0 && low < high){
				//partition
				int p;
				double pivot = V[floor((high + low)/2)].getFitness(); //get middle elm's fitness
				//std::cout << "pivot index = " << floor((high + low)/2) << std::endl;
				int i = low - 1;
				int j = high + 1;
				//std::cout << "sort() while loop..." << std::endl;
				while(true){

					do{	
						i = i+1;
					}
					while(V[i].getFitness() < pivot);	
					
					do{
						j = j-1;
					}
					while(V[j].getFitness() > pivot);

					if(i >= j){
						//std::cout << "sort: p = j" << std::endl;
						p = j; 
						break;
					}
					 
					
					//std::cout << "sort 3" << std::endl;
					Genome temp = V[i];
					//std::cout << "V[i].summary():\n";
					//V[i].summary();
					//std::cout << "temp.summary():\n";
					//temp.summary();
					//std::cout << "V[j=" << j << "].summary():";
                                        //V[j].summary();
					//std::cout << "sort 4" << std::endl;
					V[i] = V[j]; //seg fault here
					//std::cout << "sort 5" << std::endl;
					V[j] = temp;
					//std::cout << "sort: completed swap" << std::endl;
					
				}	
				//std::cout << "sort: call sort() again" << std::endl;
				sort(V, low, p);
				sort(V, p+1, high);
			}
				
		}
		
		/**
		 * @brief remove the least fit genomes from the species
		 * */
		void cull(){
			//std::cout << "inside cull() in Species.h with size = "<<  this->size()<< "\n";
			//call sort() on this species;
			//if(this->genome_vec.size() > 4){
			this->sort(this->genome_vec, 0, this->size()-1); //CHANGE: -1
			//first half of genomes and remove
			//std::cout << "[cull] after sort, just before erase\n";
			int lim = floor(this->genome_vec.size()/2);
			//std::cout << "lim=" << lim << "\n";
			//lim = 3; //TP
			this->genome_vec.erase(this->genome_vec.begin(), this->genome_vec.begin() + lim); //should erase first half
			//std::cout << "[cull]before shrinktofit...\n";
			//shrink_to_fit the genome_vec
			this->genome_vec.shrink_to_fit();
			//std::cout << "[cull]after shrinktofit...\n";
			//std::cout << "end of cull() fxn\n";
			//}
		}

	private:
		int offspring; //number of children to produce during reproduction
		double fitness; //the average (non-adjusted) fitness of this species
				
		string name;
		int num_subspecies = 0;
		int death_clock = 0; //count up each time species is sickly ( <5% size of pop and <avg pop fitness

		int last_improved = 0; //generations since this species achieved a higher fitness
		double highest_fitness = 0; //highest fitness this species has ever achieved 
};

#endif
