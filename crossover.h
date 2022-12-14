/**
 * @brief runs the NEAT algorithm on each processor
 * @author Alexander Sanfilippo
 * @date 2022-07-03
 * @TODO need crossover, and speciation functions
 *  Here is the main NEAT CPP code, which handles all
 *  parts of the program except the construction and running
 *  of the neural networks which will be done in python.  This
 *  main function will run on each processor.  
 * */


#ifndef CROSSOVER_H_7EA75F97
#define CROSSOVER_H_7EA75F97


#include <iostream>
#include <random>
#include <cstdlib>

#include "gene.h"
#include "Genome.h"
#include "Species.h"

//TEMP
#include <unistd.h>

//int Gene::innov_count{ 0 };

//using namespace std;
//int main(){

	/*Create the initial population*/

	/*MAIN LOOP*/

		/*Mutate all*/

		/*Evaluate on cart & pole*/

		/*cross-over*/




//	return 0;
//}


/**
 * @brief finds number of genomes in population.  use after cull()
 *
 * */
int calcTotalPop(std::vector <Species> *pop_ptr){
	std::vector<Species> &pop = *pop_ptr;
	int count = 0;
	for(long unsigned int i = 0; i < pop.size(); i++){
		count += pop[i].size();
		
	}		
	return count;
}

/**
 * @brief finds the population's average fitnesss
 * */
double calcAverageFitness(std::vector <Species> *pop_ptr){
	std::vector<Species> &pop = *pop_ptr;
	double avg = 0; //avg population fitness
	double pop_size = double(calcTotalPop(pop_ptr));
	//cout << "pop has " << pop_size << " genomes" << endl; 
	//double pop_size = double(pop.size());
	for(int i = 0; i < int(pop.size()); i++){ //for each species
		if(pop[i].size() > 0){ //don't count species with no genomes
			pop[i].calcFitness(); //calc species avg fitness
			double fitness = pop[i].getFitness();
			double size = pop[i].size();
			avg += fitness*size/pop_size; 
		}
		else{
			//pop_size -= 1;
		}
	}
	//avg = avg;
	return avg;
}

/**
 * @brief calculate how many offspring ALL species should make
 *
 * */

void calcOffspring(std::vector <Species> *pop_ptr, const int MAX_POP){

	/*HYPERS*/
	//if species goes 3 gens being below 5% of total pop and <half pop fitness
	//then than species is given no offspring (end of species)
	//const int MIDNIGHT = 2;
	std::vector<Species> &pop = *pop_ptr; //turn pointer into a reference.  Just a preference
	
	/*Define the number of generations a species can go without setting a new fitness high score.
 	* If this values starts at 5, then on some runs the program will kill the sole initial species
 	* before fitter members begin to appear.  To compensate, I added this block of logic to make
 	* the program more tolerant until speciation starts to occur*/
	int MAX_NO_IMPROVEMENT;
	if(pop.size() >= 8){
		MAX_NO_IMPROVEMENT = 5;
	}
	else{
		MAX_NO_IMPROVEMENT = 20;	
	
	}

	double sumOfSpeciesAdjusted = 0;
	std::vector <double> sumadj_vec;
	for(long unsigned int i = 0; i < pop.size(); i++){ //for each species
		double val = pop[i].calcSumOfAdjusted();
		sumadj_vec.push_back(val);
		//if( pop[i].size() > 1){ //don't add 1 member species into the eqn
		sumOfSpeciesAdjusted += val;
		//}
		//pop.genome_vec[i].setOffpsring( );
	}	
	//double pop_avg_fitness = calcAverageFitness(pop_ptr);
	for(long unsigned int i = 0; i < pop.size(); i++){ //for each species
			
		if(pop[i].size() > 0){
			//Elimination of Species based on gens since last improved
				
			pop[i].calcFitness(); //calc avg fitness of this species (stored in member variable
			pop[i].checkImprovement(pop[i].getFitness());//update last_improved or highest_fitness
			cout << pop[i].getLastImproved() <<" gens since " << pop[i].get_name() << " improved.\n"; //TP
			if(pop[i].getLastImproved() >= MAX_NO_IMPROVEMENT){
				pop[i].setOffspring(0); //kill the population
				std::cout << "|||||||||||||||||MAX GENS NO IMPROVEMENT HIT ON SPECIES: " << pop[i].get_name() << std::endl; 
			}
			else{ //ie, if we HAVE improved fitness in last 5 generations...
				pop[i].setOffspring( ceil( (sumadj_vec[i]/sumOfSpeciesAdjusted) * double(MAX_POP) ) );
			}
		}
		else{ //ie, if species's size is 0
			pop[i].setOffspring(0);
		}
		std::cout << "sumadj of this species = " << sumadj_vec[i] << " sumOfSpeciesAdjusted = " << sumOfSpeciesAdjusted << "\n";
		std::cout << "offspring = " << pop[i].getOffspring() << " for species "<< pop[i].get_name() <<std::endl;
	}	
}
				
				
/**
 * @breif determines compatibility of 2 genomes
 * see eqauation 1 pg 13
 * */
double compat(Genome A, Genome B){
	//find larger size of the 2 Genomes
	//std::cout << "\ninside compat()..." << std::endl;
	double D = 0; //number of disjoint genes
	double E = 0; //number of excess genes
	double N, W;
					/*Find N*/
	if(A.size() >= B.size()){ 
		N = A.size();
		//std::cout << "A is larger with size" << N << std::endl; 
	} //size of the larger of the 2 genomes
	else{ 
		N = B.size();
		//std::cout << "B is larger with size" << N << std::endl;
	}
	

					/*Find D and E genes, and W*/	
	bool a_done = false;
	bool b_done = false;
	int ia = 0;
	int ib = 0;
	double w_diff_sum = 0; //use to get W
	int comp_genes = 0; //how many compatible genes there were
	while(a_done == false || b_done == false){
		//check for match
		if(a_done == false && b_done == false){
			if(A.gene_vec[ia].getInnov() == B.gene_vec[ib].getInnov()){ //matching genes, just move on
				ia++;
				ib++;
				w_diff_sum += abs(A.gene_vec[ia].getWeight() - B.gene_vec[ib].getWeight());
				comp_genes += 1;
				//std::cout << "compate gene found" << std::endl;	
			}	
			//handle disjoints
			else{ //the innov numbers do not match
				if(A.gene_vec[ia].getInnov() > B.gene_vec[ib].getInnov()){
					D++; //we found a disjoint.
					ib++;	
					//std::cout << "disjoint gene found" << std::endl;
				}
				else{
					D++;
					ia++;
					//std::cout << "disjoint gene found" << std::endl;
				}
			}
			if(A.size() <= ia){
				a_done = true;
				//std::cout << "A hit end" << std::endl;
			}
			if(B.size() <= ib){
				b_done = true;
				//std::cout << "B hit end" << std::endl;
			}				
		}
		else{ //check for excess genes
			//std::cout << "searching for excess genes" << std::endl;
			//sleep(1); //tp
			if(a_done == true){
				E++;
				ib++;
				//std::cout << "found E on B, ib = " << ib << std::endl;
			}		
			else{
				E++;
				ia++;
				//std::cout << "found E on A" << std::endl;
			}

			if(A.size() <= ia){
				a_done = true;
			}
			if(B.size() <= ib){
				b_done = true;
			}
		
		}

	}
	
	W = w_diff_sum / comp_genes;	

	/*
	std::cout << "W = " << W << std::endl;
	std::cout << "E = " << E << std::endl;
	std::cout << "D = " << D << std::endl;
	*/
	double delta = 4*E/N + 4*D/N + 0.05*W; //1,1,0.5
	return delta;
};

/**
 * @brief combines genes of 2 parents 
 *
 * */

//void crossover(Genome *pop /*pointer to vector of Genomes*/){
	

//};
#endif

