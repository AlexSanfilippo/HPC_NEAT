/**
 * @author Alex S
 * @date 11-06-2022
 * last updated 04-08-2022
 *This file will allow us to test things as we set up the Gene, Genome, and Species classes
 * */



#include <iostream>
#include <random>

//for writing to/reading from files
#include <fstream>
#include <string>
#include <sstream>

//for calling another executable
#include <stdlib.h>

//for timing
#include <chrono>

//my headers
#include "gene.h"
#include "Genome.h"
#include "Species.h"
#include "NOV.h"
#include "crossover.h"


//REMOVE LATER
#include <typeinfo>

//MPI
#include <mpi.h>


//for python
//#include "fitness.h"


//put this function in crossover.h later
void crossover(Genome g1, Genome g2, Genome * child);
void matchMaker(Species s, int* parents, double avg_fitness);
void reproduce(std::vector <Species> *pop_ptr, const int MAX_POP, NOV *nov);
void  writeToFileRecv(ofstream *fp);


int Gene::innov_count{ 0 };

using namespace std;
int main(int argc, char **argv){	

	/*==========  Version 1.1 of Main Loop  =========*/
	//created 25-07-2022
	//updated 10-08-2022
		

	/*SIMULATION HYPERPARAMETERS*/
	const int POP_SIZE  = 150; //number of genomes in the whole population
	const int NUM_MUT = 1; //number of mutation
	const int NUM_GEN = 30; //number of generations	
	int gen_count = 0;
	NOV nov(4,1); //create Nodal Order Vector object
	
	std::vector <Species> pop; //population, a vector of species
	Species species_init = Species(0); //the initial species
	for(int i = 0; i < POP_SIZE; i++){
		Genome newGenome; //create new default genome
		species_init.genome_vec.push_back(newGenome); //add genome to the initial spec.
	}
	pop.push_back(species_init); //add species to population

	
	ofstream timefp;
	timefp.open("eval_timings.csv"); //open file for timings

	ofstream neat_time_fp;
        neat_time_fp.open("neat_timings.csv");

	ofstream fit_fp;
	fit_fp.open("pop_avg_fitness.csv");

	/*MPI EXPERIMENT*/

	int nprocs, myid;
	
	/* //old IO stuff	
	MPI_File handle;
	int buffer[3] = {7,8,9};
	*/
	MPI_Init(&argc, &argv);
	//MPI_Init(0,0);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	//std::cout << "inside "<< myid <<"/" << nprocs << " MPI processes\n"; 

	/*determine start index and quantity of genomes to write to 
 	* each file--need to do this in main loop to adjust for varying
 	* total number of genomes in the population*/
	//on all procs
	int start = 0;	
	int genomes_per_proc;
	
	std::vector <int> start_vec = {0}; //stores starting genome index for each proc
	std::vector <int> genomes_per_proc_vec; //stores ngenomes of each proc
	//on proc 0
	int ngenomes;
	if(myid == 0){
		ngenomes = calcTotalPop(&pop);
		//std::cout << "ngenomes=" << ngenomes << endl;
		//std::vector <int> genomes_per_proc; //index=proc id, stores starting index of genomes
		int small = floor(double(ngenomes)/double(nprocs));
		int big = small + 1;
		int num_big_procs = ngenomes % nprocs;
		if(num_big_procs > 0){ //proc 0 is always first to get big
			genomes_per_proc = big;
			genomes_per_proc_vec.push_back(big);
		}
		else{
			genomes_per_proc = small;
			genomes_per_proc_vec.push_back(small);
		}
		//cout << "num_big_procs=" << num_big_procs << ", small="<<small<<", big="<< big<< endl;
		for(int i = 1; i < nprocs; i++){//for each proc
			int temp_genomes_per_proc;
			if(i < num_big_procs){
				start += big;
				genomes_per_proc_vec.push_back(big);
				temp_genomes_per_proc = big;
			}
			else if( i == num_big_procs){
				start+=big;
				genomes_per_proc_vec.push_back(small);
				temp_genomes_per_proc = small;
			}
			else{
				start += small;
				genomes_per_proc_vec.push_back(small);
				temp_genomes_per_proc = small;
			}
			start_vec.push_back(start);
			MPI_Send(&start, 1, MPI_INT, i, 99, MPI_COMM_WORLD);
			MPI_Send(&temp_genomes_per_proc, 1, MPI_INT, i, 99, MPI_COMM_WORLD);
			//cout << "for proc "<<i<<" start should be " << start << endl;
		}	
	}
	if(myid != 0){
		//MPI recv to buffer "start" from rank 0
		MPI_Recv(&start, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&genomes_per_proc, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		cout << "rank " << myid << " starts at genome index " << start<<" with "<<genomes_per_proc<<" g's" << endl;
	}
	else{
		start = 0;
		cout << "rank " << myid << " starts at genome index " << start << " with "<<genomes_per_proc<<" g's" <<endl;
	}

	
	/*testing MPI-IO*/
	/*
	int inbuff[3] = {1,2,3};
	int access_mode = MPI_MODE_CREATE|MPI_MODE_WRONLY;//MPI_MODE_RDWR;
	MPI_File_open(MPI_COMM_WORLD, "io", access_mode, MPI_INFO_NULL, &handle);
	//if(myid == 0){
		if (MPI_File_write(handle, buffer, 3, MPI_INT, MPI_STATUS_IGNORE) != MPI_SUCCESS){ //write my id to file
			std::cout << "ERROR WRITING TO FILE\n";
		}
		MPI_File_close(&handle);
		std::cout << "finished write and closed the file\n";
		MPI_File_open(MPI_COMM_WORLD, "io", MPI_MODE_CREATE|MPI_MODE_RDWR, MPI_INFO_NULL, &handle);
		std::cout << "reopened the file for reading\n";
		if( MPI_File_read(handle, inbuff, 3, MPI_INT, MPI_STATUS_IGNORE) != MPI_SUCCESS){
			std::cout << "ERROR READING FROM FILE\n";
		}
		MPI_File_close(&handle);
		std::cout << "finished read and write\n";
		std::cout << "inbuff contains: " << inbuff[0]<<","<< inbuff[1]<<","<< inbuff[2] << endl;
	//}
	
	MPI_File_close(&handle);
	*/
	
	//Mutate all the genomes in the species & calc fitness
	for(int g = 0; g < NUM_GEN; g++){ //each gen

		/*[Parallel]*/
		/* //OLD
		//if the number of genomes changes, change the last procs genomes_per_proc
		int update_last_proc = 0;
		if(myid==0){
			int real_pop = calcTotalPop(&pop);
			if(ngenomes != real_pop){
				update_last_proc = real_pop - ngenomes;
				genomes_per_proc_vec[nprocs-1] += update_last_proc;
				//send signal to proc nprocs-1 to update num of genomes			
				MPI_Send(&update_last_proc, 1, MPI_INT, nprocs-1, 99, MPI_COMM_WORLD);
			}
		}
		if(myid == nprocs-1 && update_last_proc != 0){				
			MPI_Recv(&update_last_proc, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			genomes_per_proc += update_last_proc; //add genome to this proc
			//TP
			std::cout << "updated proc " <<nprocs-1 <<"'s genomes_per_proc to be " << genomes_per_proc << std::endl;
		}
		*/
		/*[Parallel] Repeat of the initial genomes/proc balancing*/
		
		
		//MPI_Barrier(MPI_COMM_WORLD);
		bool did_rebalance = false;

				
		if(myid == 0){
			int cur_ngenomes = calcTotalPop(&pop);
			if(ngenomes != cur_ngenomes){ //only rebalance if pop size has changed!
				did_rebalance = true;
			}
			
			std::cout << "~~~~~~1/2Rebalancing because old pop = " << ngenomes <<\
			" and current pop = " <<cur_ngenomes << "~~~~~~\n";
		}
		//BCAST MUST BE CALLED ON ALL PROCS AND REQS NO RECV!!!!!!!!!!!!!
		//std::cout <<"Bcasting "did_balance"\n"; 
		MPI_Bcast(&did_rebalance, 1, MPI::BOOL, 0, MPI_COMM_WORLD);
		//receive signal to rebalance on all ranks
		//MPI_Recv(&did_rebalance, 1, MPI::BOOL, 0, 9090, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Barrier(MPI_COMM_WORLD); //tp
		//std::cout << "reveived rebalance on procs\n";
		if(myid == 0 && did_rebalance == true){
			//TP
			int ngenomes_old = ngenomes;
			std::cout << "~~~~~~2/2Rebalancing because old pop = " << ngenomes_old <<\
			" and current pop = " <<ngenomes << "~~~~~~\n";
			
			int ngenomes = calcTotalPop(&pop);
			//if(ngenomes != cur_ngenomes){ //only rebalance if pop size has changed!
			//did_rebalance = true;
			//MPI_Send(&did_rebalance, 1, MPI::BOOL, i, 9090, MPI_COMM_WORLD); //tell other procs to recv new number of genomes this Gen
			//ngenomes = cur_ngenomes; //set ngenomes to current pop size	
			//clear/ reset start_vec and genomes_per_proc_vec
			start_vec.clear();
			start_vec.push_back(0);
			genomes_per_proc_vec.clear();
			//std::cout << "ngenomes=" << ngenomes << endl;
			//std::vector <int> genomes_per_proc; //index=proc id, stores starting index of genomes
			int small = floor(double(ngenomes)/double(nprocs));
			int big = small + 1;
			int num_big_procs = ngenomes % nprocs;
			if(num_big_procs > 0){ //proc 0 is always first to get big
				genomes_per_proc = big;
				genomes_per_proc_vec.push_back(big);
			}
			else{
				genomes_per_proc = small;
				genomes_per_proc_vec.push_back(small);
			}
			//cout << "num_big_procs=" << num_big_procs << ", small="<<small<<", big="<< big<< endl;
			for(int i = 1; i < nprocs; i++){//for each proc
				int temp_genomes_per_proc;
				if(i < num_big_procs){
					start += big;
					genomes_per_proc_vec.push_back(big);
					temp_genomes_per_proc = big;
				}
				else if( i == num_big_procs){
					start+=big;
					genomes_per_proc_vec.push_back(small);
					temp_genomes_per_proc = small;
				}
				else{
					start += small;
					genomes_per_proc_vec.push_back(small);
					temp_genomes_per_proc = small;
				}
				start_vec.push_back(start);
				
				MPI_Send(&start, 1, MPI_INT, i, 9091, MPI_COMM_WORLD);
				MPI_Send(&temp_genomes_per_proc, 1, MPI_INT, i, 9092, MPI_COMM_WORLD);
				//cout << "for proc "<<i<<" start should be " << start << endl;
			}	
		}
		if(myid != 0){
			if(did_rebalance == true){
				//MPI recv to buffer "start" from rank 0
				MPI_Recv(&start, 1, MPI_INT, 0, 9091, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&genomes_per_proc, 1, MPI_INT, 0, 9092, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				cout << "rank " << myid << " starts at genome index " << start<<" with "<<genomes_per_proc<<" g's" << endl;
			}
		}
		else if (myid == 0 && did_rebalance == true){
			start = 0;
			cout << "rank " << myid << " starts at genome index " << start << " with "<<genomes_per_proc<<" g's" <<endl;
		}
		
		//unsure if neccessary	
		MPI_Barrier(MPI_COMM_WORLD);

		auto start_gen_a = std::chrono::high_resolution_clock::now(); //temp fix, remove later for parallel timing
		if(myid == 0){
			start_gen_a = std::chrono::high_resolution_clock::now();
			for(long unsigned int h = 0; h < pop.size(); h++){ //each species
				for(unsigned long int i = 0; i < pop[h].genome_vec.size(); i++){ //each genome
					for(int j = 0; j < NUM_MUT; j++){ //each mutation
						pop[h].genome_vec[i].mutate(&nov);
						//std::cout << "sleeping after mutate in main()\n";
						//sleep(1);
					}
					//get fitness
					//pop[h].genome_vec[i].dummyFitness();
					//pop[h].genome_vec[i].calcAdjusted(pop[h].genome_vec.size()); //?
				}
			}
		}


		/*Write all Genomes to a File*/
		//char filename[] = "genome_data.csv";
		string filename = "genome_data";
		filename += to_string(myid) + ".csv";
		ofstream fp;
		fp.open(filename);
		int genome_count = 0; //tracks what genome we're on
		int rank = 1;
		//for each species, and each genome within, calcNodeInfo, then writeToFile
		if(myid==0){
			for(int i = 0; i < int(pop.size()); i++){//for each species 
				for(int j = 0; j < pop[i].size(); j++){//for each genome
					pop[i].genome_vec[j].calcNodeInfo(&nov);
					if(genome_count < start + genomes_per_proc){
						//std::cout << "genome " << genome_count << " written to rank " << rank<< std::endl;
						pop[i].genome_vec[j].writeToFileSend(&fp, genome_count, 0); //from sequential version
					}
					else{//if outside genomes for proc 0
						//find which rank we are writing genomes too
						bool found_proc = false;
						while(found_proc == false){
							if(genome_count < start_vec[rank] + genomes_per_proc_vec[rank]){
								//then we are in proc == rank
								/*
								std::cout << "genome_count="<<genome_count<<\
									"genome index j =" << j <<\
									",  start_vec[rank] + genomes_per_proc_vec[rank]=" <<\
									 start_vec[rank] + genomes_per_proc_vec[rank] << std::endl;
								std::cout << "sending genome to rank " << rank << std::endl;
								*/
								pop[i].genome_vec[j].writeToFileSend(&fp, genome_count, rank);
								//std::cout << "ran writeToFileSend()...\n";
								found_proc = true;
							}
							else{
								rank += 1;
							}
						}
					}
					genome_count += 1;
				}

			}
		}
		if(myid!=0){
			
			//recv weight vector from 0
			//now write this genome to local file!
			for(int i = 0; i < genomes_per_proc; i++){
				writeToFileRecv(&fp);
				
			}
		}
		
		fp.close();
			
		/*Call python script to evaluated all genome's fitnesses*/
		auto stop_gen_a = std::chrono::high_resolution_clock::now(); //stop timing before eval
		auto start = std::chrono::high_resolution_clock::now();
		
		string e_caller = "python evaluate.py " + to_string(myid);
		const char *callme = e_caller.c_str();
		//MPI_Barrier(MPI_COMM_WORLD);
		system(callme);
		MPI_Barrier(MPI_COMM_WORLD); //all ranks need to return their fitnesses first
		//const char* fitness_combiner = "cat local_fitness* > fitness.csv"
		//system(fitness_combiner);
		//system("./calls_evaluate");
		//cout << "in rank #" << myid << endl;
		//system("cat genome_data.csv");
		//std::cout <<"[MainLoop] called evaluate.py\n";
			
		auto stop = std::chrono::high_resolution_clock::now();
		auto start_gen_b = std::chrono::high_resolution_clock::now(); //start timing after eval
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		timefp << "GEN " << g << ": " << double(duration.count())/1000000 << "\n";
		//std::cout << "GEN " << g << ": " << duration.count()/1000000 << ",\n";	
		
		bool fitness_stop = false;
		/*[Parallel] Combine the fitness files into 1 big fitness file*/
		if(myid==0){
			std::cout << "rank " << myid << "combing fitnessX.csv files" << endl;
			//open another fp for writing out to file
			
			ofstream fpfitness;
			fpfitness.open("fitness.csv");

			int gcount = 0; //count which genome we're one as we go thru fitnesses
			
			for(int i = 0; i < nprocs; i++){
				string fitfile = "fitness" + to_string(i) + ".csv";
					
								
				fstream fpin; //create file pointer
				fpin.open(fitfile, ios::in); //open the file
				string line, score;
				getline(fpin, line);//reads all values and stores in "line"
				stringstream s(line);
				for(int j = 0; j <genomes_per_proc_vec[i]; j++){
					getline(s,score,',');
					//STOP CONDITION
					if(stoi(score) > 2000){
						cout << "built a perfect NN at " << gcount << endl;
						int a = 0;
						while(a <= gcount){
							for(int b = 0; b < int(pop.size()); b++){
								for(int c = 0; c < pop[b].size(); c++){
									if(a == gcount){
										cout << "FALSE-PERFECT GENOME SUMMARY:\n";
										pop[b].genome_vec[c].summary();
									}
									a += 1;
								}
							}
						}
						fitness_stop = true;
					}
					fpfitness << stoi(score) <<",";
					gcount += 1;
				}
				fpin.close();
				//read in the file data(like below)
				//write the file data to fitness.csv
			}
			fpfitness.close();
		}
		//broadcast "target fitness reached" signal to all procs from root (0)
		MPI_Bcast(&fitness_stop, 1, MPI::BOOL, 0, MPI_COMM_WORLD);
		if(myid == 0){

			/*Read the fitnesses from the file and assign to genomes' members*/
			
			std::cout << "rank 0 reading in fitness.csv into all genomes";

			fstream fpin; //create file pointer
			fpin.open("fitness.csv", ios::in); //open the file
			string line, score;
			getline(fpin, line);//reads all values and stores in "line"
			stringstream s(line);
			for(int h = 0; h < int(pop.size()); h++){//for each species
				for(int i = 0; i < int(pop[h].genome_vec.size()); i++){ //each genome
					//read fitness from the file
					getline(s,score,',');
					//do genome.setFitness(fitness)
					//cout << "read in fitness score " << score << endl;
					pop[h].genome_vec[i].setFitness(stoi(score));
					//calcAdjustedFitness	
					pop[h].genome_vec[i].calcAdjusted(pop[h].genome_vec.size()); 
				}	

			}
			fpin.close();

			/*calc species' and global average fitness scores*/
			for(int h = 0; h < int(pop.size()); h++){
				pop[h].calcFitness();
			}			


			/*=== Generational Pop. Print Out  ===*/	
			cout << "====== [pre-crossover] Generation " << gen_count << "/" << NUM_GEN-1 <<" ======\n";
			cout << "Number of Species:  " << pop.size() << "\n";
			cout << "Number of Genomes: " << calcTotalPop(&pop) << std::endl;
			double pop_avg_fitness = calcAverageFitness(&pop);
			cout << "Population Avg. Fitness (pre-cull): " << pop_avg_fitness << endl;	
			fit_fp << pop_avg_fitness << endl;	
			cout << "Species Stats:\n";
			for(int j = 0; j < (int)pop.size(); j++){
				cout << "       Species " << pop[j].get_name() << ": size: " << pop[j].size() << ", avg Fit: "\
					<< pop[j].getFitness() << "\n";
			}
	 

			
			ngenomes = calcTotalPop(&pop);// record pre-reproduce pop size for rebalancing in Next Generation
			/*reproduce (create new population from fittest of old pop)*/
			reproduce(&pop, POP_SIZE, &nov);
			cout << "[POST-CROSSOVER]Number of Genomes: " << calcTotalPop(&pop) << std::endl;
			//std::cout <<"[MainLoop] called reproduce\n";
			gen_count++;
			auto stop_gen_b = std::chrono::high_resolution_clock::now(); //stop timing before eval
			auto neat_duration_a = std::chrono::duration_cast<std::chrono::microseconds>(stop_gen_a - start_gen_a); 
			auto neat_duration_b = std::chrono::duration_cast<std::chrono::microseconds>(stop_gen_b - start_gen_b);
			neat_time_fp << "GEN " << g << ": " << double((neat_duration_a.count() + neat_duration_b.count())) / 1000000 << "\n";
			//std::cout << "nt GEN " << g << ": " << double(neat_duration_b.count()) / 1000 << "\n";
		}
		
		//if we hit a fitness score we want, stop the program
		if(fitness_stop == true){
			cout << "hit fitness_stop on rank " << myid << endl;
			MPI_Finalize();
			return 0;
		}
	}
	
	neat_time_fp.close(); //close fp for writing time of NEAT execution/Generation
	timefp.close(); //close file pointer for recording evaluation times/generation
	fit_fp.close();

	/*END OF MAIN LOOP: printing results to confirm code working*/
	/*printing results...*/
	if(myid == 0){
		cout << "\n\n\n======After Main Loop our population looks like this after "<<gen_count\
			<<" generations =======:\n";
		cout << "Number of Species:  " << pop.size() << "\n";
		cout << "Number of Genomes: " << calcTotalPop(&pop) << std::endl;
		cout << "Genomes per Species:\n";
		for(int j = 0; j < (int)pop.size(); j++){
			cout << "	Species " << pop[j].get_name() << " has " << pop[j].size() << " genomes\n";
		}
	}
	
	 //print the last two species's genomes
	/*
	cout << "Last species consists of these Genomes:\n";
	for(int i = 0; i < pop[pop.size()-1].size(); i++){ //for each genome in last species in pop
		pop[pop.size()-1].genome_vec[i].summary(); //print summary of genome
	}
				
	cout << "Second-to-Last species consists of these Genomes:\n";
	for(int i = 0; i < pop[pop.size()-2].size(); i++){ //for each genome in last species in pop
		pop[pop.size()-2].genome_vec[i].summary(); //print summary of genome
	}
	*/

	
	MPI_Finalize(); //the REAL Finalize call
	return 0;
		
}


/**
* @write to file for non-root processes
* MOVE TO TESTINGGROUNDS AND REPLACE OTHER FXN OF SAME NAME*/
void  writeToFileRecv(ofstream *fp){
	//receive weights and inputs vectors from root proc
	int node_info_size, row_lens_size; //size of inputs vector 1d
	std::vector <double> node_info_1d; //holds all connection and weight info
	std::vector <int> row_lens_1d; //holds row length of each node for a given genome
	
	MPI_Recv(&node_info_size, 1, MPI_INT, 0, 991, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//std::cout << "recv node_info_size " << node_info_size << endl;
	MPI_Recv(&row_lens_size, 1, MPI_INT, 0, 992, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//std::cout << "recv row_lens_size " << row_lens_size << endl;
	node_info_1d.reserve(node_info_size);
	row_lens_1d.reserve(row_lens_size);	
	MPI_Recv(&node_info_1d[0], node_info_size, MPI_DOUBLE, 0, 993, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//std::cout << "recv'd node_info_1d" << endl;
	MPI_Recv(&row_lens_1d[0], row_lens_size, MPI_INT, 0, 994, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//std::cout << "recv'd node_info_1d" << endl;

	//WRITE TO FILE FROM 1D:	
	int count = 0; //index for node_info_1d
	for(int i = 0; i < row_lens_size; i++){ //for each row in nodeinfo
		int row_len = row_lens_1d[i]; //get number of items in row "i"
		if(node_info_1d[count] != 0){ //if this node has at least one input
			*fp << i << ",";//write the node number first
			for(int j = 0; j < row_len; j++){ //for each item in this row
				*fp << node_info_1d[count] << ","; //write to file
				count += 1;
			}
			*fp << "\n"; //at end of row print new line 
		}
		else{
			count += row_len; //skip this whole row if it has no inputs
		}
		
	}
	*fp << "END," << "\n"; //at end of genome print "END"

}


/**
 * @brief produces offspring genome between 2 adults in a species
 * @param g1 the first parent gene
 * @param g2 the second parent gene
 * */
void crossover(Genome g1, Genome g2, Genome *child){

	//std::cout << "inside crossover()...\n";
	
	//create an empty Genome called child to store the inherited genes
	//Genome child(true);
	 
	auto urd = std::uniform_real_distribution<>(0,1); //random probability
	//auto uid = std::uniform_int_distribution<>(0,gene_vec.size() - 1); //random gene

	bool g1_done = false;
	bool g2_done = false;
	int i = 0; //index for g1's gene_vec
	int j = 0; //index for g2's gene_vec

	
	while(g1_done == false || g2_done == false){
		//check for match
		//std::cout << "b4 first if\n";
		if(g1_done == false && g2_done == false){
			//std::cout << "after first if\n";
			//std::cout << "b4 2nd if\n";
			//std:: cout << "g1.gene_vec[i].getInnov()   : "<< g1.gene_vec[i].getInnov() << std::endl;
			//std:: cout << "g2.gene_vec[j].getInnov()   : "<< g2.gene_vec[j].getInnov() << std::endl;
			if(g1.gene_vec[i].getInnov() == g2.gene_vec[j].getInnov()){ //matching genes
				//std::cout << "after 2nd if\n";
				//pick parent to inherit the gene from
				if(urd(gen) < 0.5){
					Gene gene(g1.gene_vec[i]); //cpy construct gene
					child->gene_vec.push_back(gene); //add gene to child
					
					if(child->gene_vec.back().getEnabled() == false){
						if(g1.getAdjusted() > g2.getAdjusted()){//if disabled gene has higher fitness
							if(urd(gen) > 0.75){ //25% odds to switch to enabled
								child->gene_vec.back().setEnabled(true);
							}
						}
						else{
							child->gene_vec.back().setEnabled(true);
						}
					}
				}	
				else{
					Gene gene(g2.gene_vec[j]); //cpy construct gene
                                        child->gene_vec.push_back(gene); //add gene to child
					
					if(child->gene_vec.back().getEnabled() == false){
						if(g2.getAdjusted() > g1.getAdjusted()){//if disabled gene has higher fitness
							if(urd(gen) > 0.75){ //25% odds to switch to enabled
								child->gene_vec.back().setEnabled(true);
							}
						}
						else{
							child->gene_vec.back().setEnabled(true);
						}
					}
					
				}
			
				//disabled check--OLD
				/*
				//std::cout << "b4 3rd if\n";			
				if(g1.gene_vec[i].getEnabled() == false || g2.gene_vec[j].getEnabled() == false ){
					//std::cout << "after 34d if\n";
					if(urd(gen) < 0.75){
						//set gene to disabled in child
						child->gene_vec.back().setEnabled(true);
						// *last->setEnabled(false);
					}
					else{
						//set gene to enabled in child
						child->gene_vec.back().setEnabled(true);
                                                //last.setEnabled(true);
					}
				}
				*/
				i++;
				j++;
			}	
			//handle disjoints
			else{ //the innov numbers do not match
				if(g1.gene_vec[i].getInnov() > g2.gene_vec[j].getInnov()){ //if D in g2
					//determine which parent is fitter
					if( g1.getAdjusted() <= g2.getAdjusted() ){
						//inherit the Disjoint Gene
						Gene gene ( g2.gene_vec[j] );
						child->gene_vec.push_back(gene);
					}
					
					j++;	
					//std::cout << "disjoint gene from g2 added to child" << std::endl;
				}
				else{ //if D in g1
					
					if( g1.getAdjusted() >= g2.getAdjusted() ){
						//inherit the Disjoint Gene
						Gene gene ( g1.gene_vec[i] );
						child->gene_vec.push_back(gene);
					}
					i++;
					//std::cout << "disjoint gene from g1 added to child" << std::endl;
				}
			}
			//Check for end of either Genome reached
			if(g1.size() <= i){
				g1_done = true;
				//std::cout << "g1 hit end" << std::endl;
			}
			if(g2.size() <= j){
				g2_done = true;
				//std::cout << "g2 hit end" << std::endl;
			}
						
		}
		//If one of the 2 genome's ends have been reached
		else{ //check for excess genes
			//std::cout << "searching for excess genes" << std::endl;
			//sleep(5);
			
			if(g1_done == true){ //g2 must have the E gene
				//std::cout << "adding excess gene from g2 to child...\n";
				//std::cout << "g2 has size: " << g2.size() << "and j = " << j << endl;
				if( g1.getAdjusted() <= g2.getAdjusted() ){	
					Gene temp(g2.gene_vec[j]);
					
					child->gene_vec.push_back(temp);
				}
				j++;
			}		
			else{ //g2 == done, so g1 has the E gene
				//std::cout << "adding excess gene from g1 to child...\n";
				if( g1.getAdjusted() >= g2.getAdjusted() ){
                                        Gene temp(g1.gene_vec[i]);
                                        child->gene_vec.push_back(temp);
                                }
				i++;
			}

			if(g1.size() <= i){
				g1_done = true;
				
			}
			if(g2.size() <= j){
				//std::cout << "g2.size() <= j" << endl;
				g2_done = true;
			}
	
		}
	}

	//return &child;
}

		


/**
 * @brief selects 2 parent genes within a species
 * returns the INDICES of the two parents in the species' genome_vec member
 * */
void matchMaker(Species s, int* parents, double avg_fitness){
	//int s_size = s.genome_vec.size();
	//std::cout << "[matchMaker()] with species of size " << s_size << std::endl;
	if(s.genome_vec.size() == 2){ //avoid long while loops if only 2 genomes
		//std::cout << " exactly 2 genomes in this species\n";
		parents[0] = 0;
		parents[1] = 1;
	}
	else if(s.genome_vec.size() > 2){ //make sure species has at least 2 member to breed
		//call rng
		auto uid = std::uniform_int_distribution<>(0,s.genome_vec.size() - 1); //-1 or no?
		parents[0] = uid(gen);	
		parents[1] = uid(gen);
		while(parents[0] == parents[1]){
			parents[1] = uid(gen);
		}
	}
	else if(s.genome_vec.size() < 1){
		//std::cout << " n < 2 genomes in this species" << std::endl;
		parents[0] = -1;
                parents[1] = -1;
	}
	//SEXUAL REPRODUCTION
	else if(s.genome_vec.size() == 1 && s.getFitness() > avg_fitness){
		//std::cout << "performing self-fertilization reproduction";
		parents[0] = 0;
		parents[1] = 0;
	}
	//kill unfit species with 1 genome
	else if(s.genome_vec.size() == 1 && s.getFitness() < avg_fitness){
		parents[0] = -1;
                parents[1] = -1;

	}	


}

/**
 * @brief calls all the functions related to reproduction for a given population 
 *
 * */
void reproduce(std::vector <Species> *pop_ptr, const int MAX_POP, NOV *nov){
	
	//std::cout << "in reproduce()...\n";	

	Gene ic_gene2 = Gene(0,0,0,true);
	int old_innov_count = ic_gene2.getInnovCount();
	
	//std::cout << "in reproduce() old_innov_count is " << old_innov_count << std::endl;

	std::vector<Species> &pop = *pop_ptr; //makes handling pop ptr easier

	//calc population average fitness
	double pop_avg_fitness = calcAverageFitness(pop_ptr); //crossover.h


	/*HYPER PARAMETERS*/ //add these to a struct in a world_setting.h file later on
	int dist_threshold = 3; //paper says 3
	


	//get size of population
	int pop_size = pop.size();
	calcOffspring(&pop, MAX_POP); //works for all species already, no loop req.
	//for each Species-choose a representative genome
	for(int i = 0; i < pop_size; i++){
		if(pop[i].size() > 2){ //why was this ever 4? changed to 2	
			pop[i].cull(); //ensures choose a fit rep
			pop[i].chooseRep(); //choose a random new repr genome
		}	
	}
	
	
	//NEW POPULATION-for new generation	
	std::vector <Species> new_pop;
	Species temp_spec = Species(old_innov_count);
	int sec_index = 0;
	for(long unsigned int i = 0; i < pop.size(); i++){ //for each species in old pop
		if(!(pop[i].genome_vec.empty())){ //if this old species had genomes
			Species empty_spec = Species(old_innov_count);
			new_pop.push_back(empty_spec); //add empty species to new pop
			new_pop[sec_index].rep = pop[i].rep; //deep copy?
			new_pop[sec_index].set_name( pop[i].get_name());
			new_pop[sec_index].set_num_subspecies(pop[i].get_num_subspecies());
			sec_index++;
		}
		else{
			//std::cout << "an empty species (" <<i<< ") was skipped\n";
			//sleep(2);
		}
	}
	//std::cout << "finished adding old species to new_pop, skipping empty species\n";

	int child_count = 0;

	//now make children for each species and speciate the children
	for(int i = 0; i < pop_size; i++){
		//for each offspring to make...
		for(int j = 0; j < pop[i].getOffspring(); j++){
			//matchMaker() to get two parents indices
			int parents[2] = {0,0};
			//std::cout << "calling matchMAker..."<<std::endl;
			matchMaker(pop[i], parents, pop_avg_fitness);
			//std::cout << "called matchMaker." << std::endl;
			//crossOver
			Genome child = Genome(true); //create empty child
			child_count++;
			//give child its parents' genetics
			/*std::cout << "calling crossover on species " << i << "/" << pop_size << \
				" and parents at indices " << parents[0] << ", " << parents[1]\
				<< std::endl;
			*/
			

			if(!(parents[0] == -1 || parents[1] == -1)){	
                		//std::cout << "calling crossever..."<<std::endl;
				if(parents[0] != parents[1]){
					//Genome child = Genome(true); //create empty child
					crossover(pop[i].genome_vec[parents[0]], pop[i].genome_vec[parents[1]], &child); 
				}
				//asexual reproduction v2				
				else if(parents[0] == parents[1]){ 
					/*std::cout << "attempting asexual reproduction in species " << i << "parents[0] = "<<\
						parents[0] << std::endl; */
					//just add the same identical gene to the species
					child = pop[i].genome_vec[parents[0]]; //child is CLONE of parent
					//std::cout << "cloned the child\n";
					//new_pop[i].genome_vec.push_back(child_ace); //place child in parent species
					//std::cout <<"\ntriggered asexual reproduction in species " << i << "\n";
				}

				//std::cout << "[reproduce()] after crossover()\n";
				int pop_size_dyn = new_pop.size(); //updates as new species are created;
				bool child_placed = false;
				//SPECIATE THE CHILD
				//loop through species again and loop at repr. genome to determine child's species
				for(int k = 0; k < pop_size_dyn; k++){ 
					
					double dist = compat(child/*cur child genome*/, new_pop[k].rep/*representative*/);
					if(dist <= dist_threshold){
						new_pop[k].genome_vec.push_back(child);
						child_placed = true;
						//std::cout << "Placed child " << child_count << " in species " << k << std::endl;
						k = pop_size_dyn; //stop adding the child to species
					}
				}
				
				if(child_placed == false){ //if no compatible species found for this child
					//std::cout << "adding new species..."<<std::endl;
					//then make a new species with this child as its rep
					string old_name = new_pop[i].get_name();
                                        int subspecies_count = new_pop[i].get_num_subspecies();
                                        Species new_spec = Species(child, old_innov_count, old_name, subspecies_count);  //create new species with child genome as first genome
					//above constructor handles setting the representative.
					new_pop.push_back(new_spec);
					//std::cout << "added new species."<<std::endl;
					//std::cout << "new species created for child\n";
					//std::cout << "Placed child " << child_count << " in NEW species " << std::endl;
				}
			}
			//this else block is TP only
			else {
				std::cout << "caught -1 parents at "<< i << "/" << pop_size << std::endl;
			}
		}
	}
	pop = new_pop; //TRY THIS 
	//std::cout << "reached end of reproduce()\n";
}


