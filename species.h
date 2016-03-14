/*
 Copyright 2001 The University of Texas at Austin

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
#ifndef _SPECIES_H_
#define _SPECIES_H_

#include "neat.h"
#include "organism.h"
#include "population.h"

namespace NEAT {

	class Organism;
	class Population;

	// ---------------------------------------------  
	// SPECIES CLASS:
	//   A Species is a group of similar Organisms      
	//   Reproduction takes place mostly within a
	//   single species, so that compatible organisms
	//   can mate.                                      
	// ---------------------------------------------  
	class Species {

	public:

		int id;
		int age; //The age of the Species 
		double ave_fitness; //The average fitness of the Species
		double max_fitness; //Max fitness of the Species
		double max_fitness1_ever; //The max it ever had
		double max_fitness2_ever; //The max it ever had
	        int avg_front_num; //Average of the organisms' front number
	        int num_parents;
		int expected_offspring;
		bool novel;
		bool checked;
		bool obliterate;  //Allows killing off in competitive coevolution stagnation
		std::vector<Organism*> organisms; //The organisms in the Species
		//std::vector<Organism*> reproduction_pool;  //The organisms for reproduction- NOT NEEDED 
		int age_of_last_improvement;  //If this is too long ago, the Species will goes extinct
		double average_est; //When playing real-time allows estimating average fitness

                //Size of parent population in each species. 
                //This variable is set before reproduce is called.
                //This ensures that the newly created babies are not used 
                //during reproduction (because they havent been evaluated yet)
                int poolsize;  
		bool add_Organism(Organism *o);

		Organism *first();

		bool print_to_file(std::ostream &outFile);
		bool print_to_file(std::ofstream &outFile);

		//Change the fitness of all the organisms in the species to possibly depend slightly on the age of the species
		//and then divide it by the size of the species so that the organisms in the species "share" the fitness
		void adjust_fitness();
		
                //Calculate the number of parents based on the relative species size 
                void count_parents(int total_organisms, int extra_parents);
                
                //Mark the remaining (non-parents) for elimination
                void select_parents();

                //Count the average of the organisms' front number
                void count_avg_front_num();

		double compute_average_fitness(); 

		double compute_max_fitness();

		//Counts the number of offspring expected from all its members skim is for keeping track of remaining 
		// fractional parts of offspring and distributing them among species
		double count_offspring(double skim);

		//Compute generations since last improvement
		int last_improved() {
			return age-age_of_last_improvement;
		}

		//Remove an organism from Species
		bool remove_org(Organism *org);

		double size() {
			return organisms.size();
		}

		Organism *get_champ();

		//Perform mating and mutation to form next generation
		bool reproduce(int generation, Population *pop,std::vector<Species*> &sorted_species);

		//Perform mating and mutation to form next generation (Aditya-NSGA-2)
		bool reproduce_multiobj(int generation, Population *pop);
		
                //Binary tournament selection
                Organism* binary_tournament_select(int size);
                
                //Compare two organisms based NSGA-2 definition
                bool org1_wins_org2(Organism* org1, Organism* org2);
                
                // *** Real-time methods *** 

		//Place organisms in this species in order by their orig fitness 1
		bool rank_orig_fitness1();
		
                //Place organisms in this species in order by their orig fitness 2
		bool rank_orig_fitness2();

                //Place organisms in this species in order by their front number and crowding distance 
		bool rank_front_num_crowd_dist();

		Species(int i);

		//Allows the creation of a Species that won't age (a novel one)
		//This protects new Species from aging inside their first generation
		Species(int i,bool n);

		~Species();

	};

	// This is used for list sorting of Species by fitness of best organism highest fitness first 
	bool order_species(Species *x, Species *y);
	
	// This is used for list sorting of Species by front_num of best organism.. smallest front_num first 
        bool order_species_by_front_num_crowd_dist(Species *x, Species *y);

	bool order_new_species(Species *x, Species *y);

}

#endif
