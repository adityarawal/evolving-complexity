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
#include "species.h"
#include "organism.h"
#include <cmath>
#include <iostream>
using namespace NEAT;

Species::Species(int i) {
	id=i;
	age=1;
	ave_fitness=0.0;
	expected_offspring=0;
	novel=false;
	age_of_last_improvement=0;
	max_fitness=0;
	max_fitness1_ever=0;
	max_fitness2_ever=0;
	avg_front_num = 0;
        obliterate=false;
        poolsize = 0;

	average_est=0;
}

Species::Species(int i,bool n) {
	id=i;
	age=1;
	ave_fitness=0.0;
	expected_offspring=0;
	novel=n;
	age_of_last_improvement=0;
	max_fitness=0;
	max_fitness1_ever=0;
	max_fitness2_ever=0;
	avg_front_num = 0;
	obliterate=false;
        poolsize = 0;

	average_est=0;
}


Species::~Species() {

	std::vector<Organism*>::iterator curorg;

	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		delete (*curorg);
	}

}

bool Species::rank_orig_fitness1() {
	//organisms.qsort(order_orgs);
        std::sort(organisms.begin(), organisms.end(), order_orgs_by_orig_fitness1);
	return true;
}

bool Species::rank_orig_fitness2() {
	//organisms.qsort(order_orgs);
        std::sort(organisms.begin(), organisms.end(), order_orgs_by_orig_fitness2);
	return true;
}

bool Species::rank_front_num_crowd_dist() {
        std::sort(organisms.begin(), organisms.end(), order_orgs_by_front_num_crowd_dist);
	return true;
}

bool Species::add_Organism(Organism *o){
	organisms.push_back(o);
	return true;
}

Organism *Species::get_champ() {
	double champ_fitness=-1.0;
	Organism *thechamp;
	std::vector<Organism*>::iterator curorg;

	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		//TODO: Remove DEBUG code
		//cout<<"searching for champ...looking at org "<<(*curorg)->gnome->genome_id<<" fitness: "<<(*curorg)->fitness<<endl;
		if (((*curorg)->fitness1)>champ_fitness) {
			thechamp=(*curorg);
			champ_fitness=thechamp->fitness1;
		}
	}

	//cout<<"returning champ #"<<thechamp->gnome->genome_id<<endl;

	return thechamp;

}

bool Species::remove_org(Organism *org) {
	std::vector<Organism*>::iterator curorg;

	curorg=organisms.begin();
	while((curorg!=organisms.end())&&
		((*curorg)!=org))
		++curorg;

	if (curorg==organisms.end()) {
		//cout<<"ALERT: Attempt to remove nonexistent Organism from Species"<<endl;
		return false;
	}
	else {
		organisms.erase(curorg);
		return true;
	}

}

Organism *Species::first() {
	return *(organisms.begin());
}
/*
bool Species::print_to_file(std::ostream &outFile) {
	std::vector<Organism*>::iterator curorg;

	//Print a comment on the Species info
	//outFile<<endl<<"/* Species #"<<id<<" : (Size "<<organisms.size()<<") (AF "<<ave_fitness<<") (Age "<<age<<")  *///"<<endl<<endl;
	//char tempbuf[1024];
	//sprintf(tempbuf, sizeof(tempbuf), "/* Species #%d : (Size %d) (AF %f) (Age %d)  */\n\n", id, organisms.size(), average_est, age);
	//sprintf(tempbuf, sizeof(tempbuf), "/* Species #%d : (Size %d) (AF %f) (Age %d)  */\n\n", id, organisms.size(), ave_fitness, age);
	//outFile.write(strlen(tempbuf), tempbuf);

	//Show user what's going on
	//cout<<endl<<"/* Species #"<<id<<" : (Size "<<organisms.size()<<") (AF "<<ave_fitness<<") (Age "<<age<<")  */"<<endl;

	//Print all the Organisms' Genomes to the outFile
	//for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {

		//Put the fitness for each organism in a comment
		//outFile<<endl<<"/* Organism #"<<((*curorg)->gnome)->genome_id<<" Fitness: "<<(*curorg)->fitness<<" Error: "<<(*curorg)->error<<" */"<<endl;

	//	char tempbuf2[1024];
	//	sprintf(tempbuf2, sizeof(tempbuf2), "/* Organism #%d Fitness: %f Error: %f */\n", ((*curorg)->gnome)->genome_id, (*curorg)->fitness, (*curorg)->error);
	//	outFile.write(strlen(tempbuf2), tempbuf2);

		//If it is a winner, mark it in a comment
	//	if ((*curorg)->winner) {
	//		char tempbuf3[1024];
	//		sprintf(tempbuf3, sizeof(tempbuf3), "/* ##------$ WINNER %d SPECIES #%d $------## */\n", ((*curorg)->gnome)->genome_id, id);
			//outFile<<"/* ##------$ WINNER "<<((*curorg)->gnome)->genome_id<<" SPECIES #"<<id<<" $------## */"<<endl;
	//	}

	//	((*curorg)->gnome)->print_to_file(outFile);
		//We can confirm by writing the genome #'s to the screen
		//cout<<((*curorg)->gnome)->genome_id<<endl;
	//}

	//return true;

//}*/

//Print Species to a file outFile
bool Species::print_to_file(std::ofstream &outFile) {
  std::vector<Organism*>::iterator curorg;

  //Print a comment on the Species info
  outFile<<std::endl<<"/* Species #"<<id<<" : (Size "<<organisms.size()<<") (AF "<<ave_fitness<<") (Age "<<age<<")  */"<<std::endl<<std::endl;

  //Show user what's going on
  std::cout<<std::endl<<"/* Species #"<<id<<" : (Size "<<organisms.size()<<") (AF "<<ave_fitness<<") (Age "<<age<<")  */"<<std::endl;

  //Print all the Organisms' Genomes to the outFile
  for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {

    //Put the fitness for each organism in a comment
    outFile<<std::endl<<"/* Organism #"<<((*curorg)->gnome)->genome_id<<" Fitness1: "<<(*curorg)->fitness1<<" Fitness2: "<<(*curorg)->fitness2<<" Error: "<<(*curorg)->error<<" */"<<std::endl;

    //If it is a winner, mark it in a comment
    if ((*curorg)->winner) outFile<<"/* ##------$ WINNER "<<((*curorg)->gnome)->genome_id<<" SPECIES #"<<id<<" $------## */"<<std::endl;

    ((*curorg)->gnome)->print_to_file(outFile);
    //We can confirm by writing the genome #'s to the screen
    //std::cout<<((*curorg)->gnome)->genome_id<<std::endl;
  }

  return true;

}


bool Species::print_to_file(std::ostream &outFile) {
	std::vector<Organism*>::iterator curorg;

	//Print a comment on the Species info
	//outFile<<std::endl<<"/* Species #"<<id<<" : (Size "<<organisms.size()<<") (AF "<<ave_fitness<<") (Age "<<age<<")  */"<<std::endl<<std::endl;
	char tempbuf[1024];
	sprintf(tempbuf,"/* Species #%d : (Size %d) (AF %f) (Age %d)  */\n\n", id, organisms.size(), average_est, age);
	//sprintf(tempbuf, "/* Species #%d : (Size %d) (AF %f) (Age %d)  */\n\n", id, organisms.size(), ave_fitness, age);
	outFile << tempbuf;

	//Show user what's going on
	//std::cout<<std::endl<<"/* Species #"<<id<<" : (Size "<<organisms.size()<<") (AF "<<ave_fitness<<") (Age "<<age<<")  */"<<std::endl;

	//Print all the Organisms' Genomes to the outFile
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {

		//Put the fitness for each organism in a comment
		//outFile<<std::endl<<"/* Organism #"<<((*curorg)->gnome)->genome_id<<" Fitness: "<<(*curorg)->fitness<<" Error: "<<(*curorg)->error<<" */"<<std::endl;
		char tempbuf2[1024];
		sprintf(tempbuf2, "/* Organism #%d Fitness1: %d Fitness2: %f Time: %d */\n", ((*curorg)->gnome)->genome_id, (*curorg)->fitness1, (*curorg)->fitness2, (*curorg)->time_alive);
		outFile << tempbuf2;

		//If it is a winner, mark it in a comment
		if ((*curorg)->winner) {
			char tempbuf3[1024];
			sprintf(tempbuf3, "/* ##------$ WINNER %d SPECIES #%d $------## */\n", ((*curorg)->gnome)->genome_id, id);
			//outFile<<"/* ##------$ WINNER "<<((*curorg)->gnome)->genome_id<<" SPECIES #"<<id<<" $------## */"<<std::endl;
		}

		((*curorg)->gnome)->print_to_file(outFile);
		//We can confirm by writing the genome #'s to the screen
		//std::cout<<((*curorg)->gnome)->genome_id<<std::endl;
	}
	char tempbuf4[1024];
	sprintf(tempbuf4, "\n\n");
	outFile << tempbuf4;

	return true;

}


//Prints the champions of each species to files    
//starting with directory_prefix
//The file name are as follows: [prefix]g[generation_num]cs[species_num]
//Thus, they can be indexed by generation or species
//bool Population::print_species_champs_tofiles(char *directory_prefix, int generation) {
//
//ostrstream *fnamebuf; //File for output
//std::vector<Species*>::iterator curspecies;
//Organism *champ;
//int pause;
//
//std::cout<<generation<<std::endl;
//std::cout<<"Printing species champs to file"<<std::endl;
////cin>>pause;
//
////Step through the Species and print their champs to files
//for(curspecies=species.begin();curspecies!=species.end();++curspecies) {
//
//std::cout<<"Printing species "<<(*curspecies)->id<<" champ to file"<<std::endl;
//
////cin>>pause;
//
////Get the champ of this species
//champ=(*curspecies)->get_champ();
//
////Revise the file name
//fnamebuf=new ostrstream();
//(*fnamebuf)<<directory_prefix<<"g"<<generation<<"cs"<<(*curspecies)->id<<ends;  //needs end marker
//
////Print to file using organism printing (includes comments)
//champ->print_to_file(fnamebuf->str());
//
////Reset the name
//fnamebuf->clear();
//delete fnamebuf;
//}
//return true;
//}

void Species::adjust_fitness() {
	std::vector<Organism*>::iterator curorg;

	int count;

	int age_debt; 

	//std::cout<<"Species "<<id<<" last improved "<<(age-age_of_last_improvement)<<" steps ago when it moved up to "<<max_fitness_ever<<std::endl;

	age_debt=(age-age_of_last_improvement+1)-NEAT::dropoff_age;

	if (age_debt==0) age_debt=1;

	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {

		//Remember the original fitness before it gets modified
		(*curorg)->orig_fitness1=(*curorg)->fitness1;
		(*curorg)->orig_fitness2=(*curorg)->fitness2;

		//Make fitness decrease after a stagnation point dropoff_age
		//Added an if to keep species pristine until the dropoff point
		//obliterate is now used to make sure that num_parents for this 
                //this species = 0 i.e. the species does not reproduce and dies 
		//Low fitness reduces the number of expected_offsprings
		if ((age_debt>=1)) {
                        obliterate = true; //New way: obliterate marks the elimination of the species

			//Possible graded dropoff
			//((*curorg)->fitness)=((*curorg)->fitness)*(-atan(age_debt));

			//Extreme penalty for a long period of stagnation (divide fitness by 100)
			((*curorg)->fitness1)=((*curorg)->fitness1)*0.01;
			((*curorg)->fitness2)=((*curorg)->fitness2)*0.01;
			//std::cout<<"OBLITERATE Species "<<id<<" of age "<<age<<std::endl;
			//std::cout<<"dropped fitness to "<<((*curorg)->fitness)<<std::endl;
		}

		//Give a fitness boost up to some young age (niching)
		//The age_significance parameter is a system parameter
		//  if it is 1, then young species get no fitness boost
		if (age<=10) {
                        ((*curorg)->fitness1)=((*curorg)->fitness1)*NEAT::age_significance; 
                        ((*curorg)->fitness2)=((*curorg)->fitness2)*NEAT::age_significance;
                } 

		//Do not allow negative fitness
		if (((*curorg)->fitness1)<0.0) (*curorg)->fitness1=0.0001; 
		if (((*curorg)->fitness2)<0.0) (*curorg)->fitness2=0.0001; 

		//Share fitness with the species
		(*curorg)->fitness1=((*curorg)->fitness1)/(organisms.size());
		(*curorg)->fitness2=((*curorg)->fitness2)/(organisms.size());

	}

	//Sort the population by fitness1
	std::sort(organisms.begin(), organisms.end(), order_orgs_by_orig_fitness1);

	//Update age_of_last_improvement for Fitness 1 here
	if (((*(organisms.begin()))->orig_fitness1)> 
	    max_fitness1_ever) {
	  age_of_last_improvement=age;
	  max_fitness1_ever=((*(organisms.begin()))->orig_fitness1);
	}
	
        //Sort the population by fitness2
	std::sort(organisms.begin(), organisms.end(), order_orgs_by_orig_fitness2);

	//Update age_of_last_improvement for Fitness 2 here
	if (((*(organisms.begin()))->orig_fitness2)> 
	    max_fitness2_ever) {
	  age_of_last_improvement=age;
	  max_fitness2_ever=((*(organisms.begin()))->orig_fitness2);
	}

}

void Species::count_avg_front_num() {
        
	std::vector<Organism*>::iterator curorg;
        avg_front_num = 0;
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
                avg_front_num += (*curorg)->front_num;
        }
        avg_front_num = avg_front_num/organisms.size();
}

void Species::count_parents(int total_organisms, int extra_parents) {

        double relative_size;
        //Old: Adding 1.0 ensures that at least one will survive
        //Until obliterate is true (NEAT::dropoff_age), the species keeps surviving with at least one organism
        
        //New: Removing +1.0 (Minimum one parent is guranteed to each species previously)
        //Number of parents in a species depends on the relative species size
        if (!obliterate) {
                relative_size = ((double)organisms.size()/(double)total_organisms);
                num_parents += (int) floor(relative_size*extra_parents); 
        }
        else {
                num_parents = 0; 
        }
	
}

void Species::select_parents() {

	std::vector<Organism*>::iterator curorg;
        int count;
      
        //If species is marked for death, set eliminate flag for each organism and return 
        if (obliterate) {
                for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
                        (*curorg)->eliminate = true;
                }
                return;
        }

	//Mark for death those who are ranked too low to be parents 
        //Species has been already sorted by front num previously
	curorg=organisms.begin();
	//(*curorg)->champion=true;  //Mark the champ as such
	for(count=1;count<=num_parents;count++) {
	  if (curorg!=organisms.end())
	    ++curorg;
	}
	while(curorg!=organisms.end()) {
	  (*curorg)->eliminate=true;  //Mark for elimination
	  std::cout<<"marked org # "<<(*curorg)->gnome->genome_id<<" fitness = "<<(*curorg)->fitness1<<std::endl;
	  ++curorg;
	}             
        
}

double Species::compute_average_fitness() {
	std::vector<Organism*>::iterator curorg;

	double total=0.0;

	//int pause; //DEBUG: Remove

	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		total+=(*curorg)->fitness1;
		//std::cout<<"new total "<<total<<std::endl; //DEBUG: Remove
	}

	ave_fitness=total/(organisms.size());

	//DEBUG: Remove
	//std::cout<<"average of "<<(organisms.size())<<" organisms: "<<ave_fitness<<std::endl;
	//cin>>pause;

	return ave_fitness;

}

double Species::compute_max_fitness() {
	double max=0.0;
	std::vector<Organism*>::iterator curorg;

	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		if (((*curorg)->fitness1)>max)
			max=(*curorg)->fitness1;
	}

	max_fitness=max;

	return max;
}

double Species::count_offspring(double skim) {
	std::vector<Organism*>::iterator curorg;
	int e_o_intpart;  //The floor of an organism's expected offspring
	double e_o_fracpart; //Expected offspring fractional part
	double skim_intpart;  //The whole offspring in the skim

	expected_offspring=0;

	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		e_o_intpart=(int) floor((*curorg)->expected_offspring);
		e_o_fracpart=fmod((*curorg)->expected_offspring,1.0);

		expected_offspring+=e_o_intpart;

		//Skim off the fractional offspring
		skim+=e_o_fracpart;

		//NOTE:  Some precision is lost by computer
		//       Must be remedied later
		if (skim>1.0) {
			skim_intpart=floor(skim);
			expected_offspring+=(int) skim_intpart;
			skim-=skim_intpart;
		}
	}

	return skim;

}

bool Species::org1_wins_org2(Organism* org1, Organism* org2) {
        if ((org1->front_num < 0) ||(org2->front_num < 0))  {//For verification purpose
                std::cout<<"ERROR: Incorrect front number"<<org1->front_num<<" "<< org2->front_num<<std::endl;
                exit(0);
        }
        if (org1->front_num < org2->front_num) {//Lower numbered fronts are more dominant
                return true;
        }
        else if (org1->front_num > org2->front_num) {
                return false;
        }
        else if (org1->front_num == org2->front_num) {
                if (org1->crowd_dist > org2->crowd_dist) {
                        return true;
                }
                else if (org1->crowd_dist < org2->crowd_dist) {
                        return false;
                }
                else if (org1->crowd_dist == org2->crowd_dist){//Doesn't matter whether we return true or false
                        return true;
                }
        }
                
}

Organism* Species::binary_tournament_select(int size) {//Passing size as argument is important because species size increases as new babies
                                                       //are added. Need to make sure that we pass the original species size
	
        std::vector<Organism*>::iterator curorg;
        Organism* org1; 
        Organism* org2;
        int orgnum;
        int orgcount;

        //RANDOM PARENT CHOOSER
	orgnum=rand()%size; //Random integer between 0 and size-1
        curorg=organisms.begin();
	for(orgcount=0;orgcount<orgnum;orgcount++){
		++curorg;                       
        }
        org1 = (*curorg);//organisms[orgnum];
        //RANDOM PARENT CHOOSER
	orgnum=rand()%size; //Random integer between 0 and size-1
	curorg=organisms.begin();
	for(orgcount=0;orgcount<orgnum;orgcount++){
		++curorg;                       
        }
        org2 = (*curorg);//organisms[orgnum];

        if (org1_wins_org2(org1, org2)) {
                return org1;
        }
        else {
                return org2; 
        }
}
void check_lstm_genes(Genome *last_genome, Genome *new_genome, char *s){//Added for verification purposes only

	std::vector<Gene*>::iterator p1gene;
	std::vector<Gene*>::iterator p2gene;
        for (p1gene=(new_genome->genes).begin(); p1gene!=(new_genome->genes).end(); p1gene++) {
                if (!((*p1gene)->gate_type==(int)(*p1gene)->mutation_num)) {
                std::cout<<s<<" gate_type does not match mutation_num: "<<(*p1gene)->innovation_num<<" "<<(*p1gene)->gate_type<<" "<<(*p1gene)->mutation_num<<std::endl;
                                                std::ofstream outFile("last_genome.txt",std::ios::out);
	                                        last_genome->print_to_file(outFile);
                                                outFile.close();
                                                std::ofstream outFile2("new_genome.txt",std::ios::out);
	                                        new_genome->print_to_file(outFile2);
                                                outFile2.close();
                exit(0);
                }
        }
        for (p1gene=(last_genome->genes).begin(); p1gene!=(last_genome->genes).end(); p1gene++) {
                for (p2gene=(new_genome->genes).begin(); p2gene!=(new_genome->genes).end(); p2gene++) {
                        if ((*p1gene)->innovation_num == (*p2gene)->innovation_num) {
                                if (!((*p1gene)->gate_type == (*p2gene)->gate_type)) {
                                                std::cout<<s<<" Gate_type mismatch ERROR  ";
                                                std::cout<<(*p1gene)->innovation_num<<" "<<(*p1gene)->gate_type<<" "<< (*p2gene)->gate_type <<std::endl;
                                                std::cout<<(*p1gene)->innovation_num<<" "<<(*p1gene)->mutation_num<<" "<< (*p2gene)->mutation_num <<std::endl;
                                                std::ofstream outFile("last_genome.txt",std::ios::out);
	                                        last_genome->print_to_file(outFile);
                                                outFile.close();
                                                std::ofstream outFile2("new_genome.txt",std::ios::out);
	                                        new_genome->print_to_file(outFile2);
                                                outFile2.close();
                                        exit(0);
                                }
                                if (!((*p1gene)->mutation_num == (*p2gene)->mutation_num)) {
                                                std::cout<<s<<" Mutation_num mismatch ERROR  ";
                                                std::cout<<(*p1gene)->innovation_num<<" "<<(*p1gene)->mutation_num<<" "<< (*p2gene)->mutation_num <<std::endl;
                                                std::ofstream outFile("last_genome.txt",std::ios::out);
	                                        last_genome->print_to_file(outFile);
                                                outFile.close();
                                                std::ofstream outFile2("new_genome.txt",std::ios::out);
	                                        new_genome->print_to_file(outFile2);
                                                outFile2.close();
                                        exit(0);
                                }
                        }
                }
        }
}
bool Species::reproduce_multiobj(int generation, Population *pop) {//(Aditya - NSGA2)
	int count;
	std::vector<Organism*>::iterator curorg;

	//int poolsize;  //The number of Organisms in the old generation

	Organism *mom; //Parent Organisms
	Organism *dad;
	Organism *baby;  //The new Organism

	Genome *new_genome;  //For holding baby's genes

	std::vector<Species*>::iterator curspecies;  //For adding baby

	double randmult;

	Network *net_analogue;  //For adding link to test for recurrency
	int pause;

	bool outside;

	bool found;  //When a Species is found

	bool champ_done=false; //Flag the preservation of the champion  

	int giveup; //For giving up finding a mate outside the species

	bool mut_struct_baby;
	bool mate_baby;

	//The weight mutation power is species specific depending on its age
	double mut_power=NEAT::weight_mut_power;

        if (organisms.size()!= NEAT::pop_size/2) {//For verification purposes
                std::cout<<"ERROR- mismatch in organism size and NEAT::pop_size/2"<<std::endl;
                exit(0);
        }
        expected_offspring = organisms.size();//Number of new organisms that will be generated here
        poolsize = organisms.size();//Size of the parent population used for binary tournament selection

        for (count= 0; count <expected_offspring; count++) { //Create NEAT::pop_size/2 new offsprings
		
		mut_struct_baby=false;
		mate_baby=false;
		outside=false;


                //Mutate only 
		if (randfloat()<NEAT::mutate_only_prob) {

			//Choose the parent
                        mom = binary_tournament_select(poolsize);

			new_genome=(mom->gnome)->duplicate(count);

			//Do the mutation depending on probabilities of 
			//various mutations

			if (randfloat()<NEAT::mutate_add_node_prob) {
			        Genome *last_genome=(new_genome)->duplicate(1);
				//std::cout<<"mutate add node"<<std::endl;
				new_genome->mutate_add_node(pop->innovations,pop->cur_node_id,pop->cur_innov_num);
				mut_struct_baby=true;
                                char *s = "mutate_add_node";
                                check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
			}
			else if (randfloat()<NEAT::mutate_add_lstm_node_prob) {//Adding LSTM node is less preferred than regular node
			        Genome *last_genome=(new_genome)->duplicate(1);
                                new_genome->mutate_add_lstm_node(pop->innovations,pop->cur_node_id,pop->cur_innov_num);
				mut_struct_baby=true;
                                char *s = "mutate_add_lstm_node";
                                check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
                        }
			else if (randfloat()<NEAT::mutate_add_link_prob) {
				//std::cout<<"mutate add link"<<std::endl;
				net_analogue=new_genome->genesis(generation);
			        Genome *last_genome=(new_genome)->duplicate(1);
				new_genome->mutate_add_link(pop->innovations,pop->cur_innov_num,NEAT::newlink_tries);
                                char *s = "mutate_add_link";
                                check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
				delete net_analogue;
				mut_struct_baby=true;
			}
			//NOTE:  A link CANNOT be added directly after a node was added because the phenotype
			//       will not be appropriately altered to reflect the change
			else {
				//If we didn't do a structural mutation, we do the other kinds

				if (randfloat()<NEAT::mutate_random_trait_prob) {
					//std::cout<<"mutate random trait"<<std::endl;
					new_genome->mutate_random_trait();
				}
				if (randfloat()<NEAT::mutate_link_trait_prob) {
					//std::cout<<"mutate_link_trait"<<std::endl;
					new_genome->mutate_link_trait(1);
				}
				if (randfloat()<NEAT::mutate_node_trait_prob) {
					//std::cout<<"mutate_node_trait"<<std::endl;
					new_genome->mutate_node_trait(1);
				}
				if (randfloat()<NEAT::mutate_link_weights_prob) {
					//std::cout<<"mutate_link_weights"<<std::endl;
			                Genome *last_genome=(new_genome)->duplicate(1);
					new_genome->mutate_link_weights(mut_power,1.0,GAUSSIAN);
                                        char *s = "mutate_link_weights";
                                        check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
				}
				if (randfloat()<NEAT::mutate_toggle_enable_prob) {
					//std::cout<<"mutate toggle enable"<<std::endl;
			                Genome *last_genome=(new_genome)->duplicate(1);
					new_genome->mutate_toggle_enable(1);
                                        char *s = "mutate_add_toggle_enable";
                                        check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
				}
				if (randfloat()<NEAT::mutate_gene_reenable_prob) {
					//std::cout<<"mutate gene reenable"<<std::endl;
			                Genome *last_genome=(new_genome)->duplicate(1);
					new_genome->mutate_gene_reenable();
                                        char *s = "mutate_add_gene_reenable";
                                        check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
				}
			}

			baby=new Organism(0.0, 0.0, new_genome,generation);
                }

                //Crossover and Mutate
                else {
			//Choose parents
                        mom = binary_tournament_select(poolsize);
                        dad = binary_tournament_select(poolsize);

			//Perform mating based on probabilities of differrent mating types
			if (randfloat()<NEAT::mate_multipoint_prob) { 
			        Genome *last_genome=(mom->gnome)->duplicate(1);
				new_genome=(mom->gnome)->mate_multipoint(dad->gnome,count,mom->front_num,dad->front_num, mom->crowd_dist,dad->crowd_dist,outside);
                                char *s = "mutate_mate_multipoint";
                                check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
			}
			else if (randfloat()<(NEAT::mate_multipoint_avg_prob/(NEAT::mate_multipoint_avg_prob+NEAT::mate_singlepoint_prob))) {
			        Genome *last_genome=(mom->gnome)->duplicate(1);
				new_genome=(mom->gnome)->mate_multipoint_avg(dad->gnome,count,mom->front_num,dad->front_num, mom->crowd_dist,dad->crowd_dist,outside);
                                char *s = "mutate_mate_multipoint_avg";
                                check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
			}
			else {
			        Genome *last_genome=(mom->gnome)->duplicate(1);
				new_genome=(mom->gnome)->mate_singlepoint(dad->gnome,count);
                                char *s = "mate_singlepoint";
                                check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
			}

			mate_baby=true;

			//Determine whether to mutate the baby's Genome
			//This is done randomly or if the mom and dad are the same organism
			if ((randfloat()>NEAT::mate_only_prob)||
				((dad->gnome)->genome_id==(mom->gnome)->genome_id)||
				(((dad->gnome)->compatibility(mom->gnome))==0.0))
			{
				//Do the mutation depending on probabilities of 
				//various mutations
				if (randfloat()<NEAT::mutate_add_node_prob) {
			                Genome *last_genome=(new_genome)->duplicate(1);
					new_genome->mutate_add_node(pop->innovations,pop->cur_node_id,pop->cur_innov_num);
					//  std::cout<<"mutate_add_node: "<<new_genome<<std::endl;
					mut_struct_baby=true;
                                        char *s = "mutate_add_node";
                                        check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
				}
		        	else if (randfloat()<NEAT::mutate_add_lstm_node_prob) {//Adding LSTM node is less preferred than regular node
			                Genome *last_genome=(new_genome)->duplicate(1);
                                        new_genome->mutate_add_lstm_node(pop->innovations,pop->cur_node_id,pop->cur_innov_num);
		        		mut_struct_baby=true;
                                        char *s = "mutate_add_lstm_node";
                                        check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
                                }
				else if (randfloat()<NEAT::mutate_add_link_prob) {
					net_analogue=new_genome->genesis(generation);
			                Genome *last_genome=(new_genome)->duplicate(1);
					new_genome->mutate_add_link(pop->innovations,pop->cur_innov_num,NEAT::newlink_tries);
                                        char *s = "mutate_add_link";
                                        check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
					delete net_analogue;
					//std::cout<<"mutate_add_link: "<<new_genome<<std::endl;
					mut_struct_baby=true;
				}
				else {
					//Only do other mutations when not doing sturctural mutations

					if (randfloat()<NEAT::mutate_random_trait_prob) {
						new_genome->mutate_random_trait();
						//std::cout<<"..mutate random trait: "<<new_genome<<std::endl;
					}
					if (randfloat()<NEAT::mutate_link_trait_prob) {
						new_genome->mutate_link_trait(1);
						//std::cout<<"..mutate link trait: "<<new_genome<<std::endl;
					}
					if (randfloat()<NEAT::mutate_node_trait_prob) {
						new_genome->mutate_node_trait(1);
						//std::cout<<"mutate_node_trait: "<<new_genome<<std::endl;
					}
					if (randfloat()<NEAT::mutate_link_weights_prob) {
			                        Genome *last_genome=(new_genome)->duplicate(1);
						new_genome->mutate_link_weights(mut_power,1.0,GAUSSIAN);
						//std::cout<<"mutate_link_weights: "<<new_genome<<std::endl;
                                                char *s = "mutate_link_weights";
                                                check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
					}
					if (randfloat()<NEAT::mutate_toggle_enable_prob) {
			                        Genome *last_genome=(new_genome)->duplicate(1);
						new_genome->mutate_toggle_enable(1);
						//std::cout<<"mutate_toggle_enable: "<<new_genome<<std::endl;
                                                char *s = "mutate_toggle_enable";
                                                check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
					}
					if (randfloat()<NEAT::mutate_gene_reenable_prob) {
			                        Genome *last_genome=(new_genome)->duplicate(1);
						new_genome->mutate_gene_reenable(); 
						//std::cout<<"mutate_gene_reenable: "<<new_genome<<std::endl;
                                                char *s = "mutate_gene_reenable_prob";
                                                check_lstm_genes(last_genome, new_genome, s);
                                delete last_genome;
					}
				}

				//Create the baby
				baby=new Organism(0.0, 0.0, new_genome,generation);

			}
			else {
				//Create the baby without mutating first
				baby=new Organism(0.0, 0.0, new_genome,generation);
			}

		}

		//Add the baby to its proper Species
		//If it doesn't fit a Species, create a new one
		baby->mut_struct_baby=mut_struct_baby;
		baby->mate_baby=mate_baby;

		curspecies=(pop->species).begin();//Only one species in NSGA-2
		(*curspecies)->add_Organism(baby);
		baby->species=(*curspecies);  //Point organism to its species

        }   

	return true;
}

bool Species::reproduce(int generation, Population *pop,std::vector<Species*> &sorted_species) {
	int count;
	std::vector<Organism*>::iterator curorg;

	int orgnum;  //Random variable
	int orgcount;
	Organism *mom; //Parent Organisms
	Organism *dad;
	Organism *baby;  //The new Organism

	Genome *new_genome;  //For holding baby's genes

	std::vector<Species*>::iterator curspecies;  //For adding baby
	Species *newspecies; //For babies in new Species
	Organism *comporg;  //For Species determination through comparison

	Species *randspecies;  //For mating outside the Species
	double randmult;
	int randspeciesnum;
	int spcount;  
	std::vector<Species*>::iterator cursp;

	Network *net_analogue;  //For adding link to test for recurrency
	int pause;

	bool outside;

	bool found;  //When a Species is found

	bool champ_done=true; //Do not need to duplicate the champ in elitist NSGA-2. Flag the preservation of the champion  

	Organism *thechamp;

	int giveup; //For giving up finding a mate outside the species

	bool mut_struct_baby;
	bool mate_baby;

	//The weight mutation power is species specific depending on its age
	double mut_power=NEAT::weight_mut_power;

	//Roulette wheel variables
	double total_fitness=0.0;
	double marble;  //The marble will have a number between 0 and total_fitness
	double spin;  //0Fitness total while the wheel is spinning

	//Compute total fitness of species for a roulette wheel
	//Note: You don't get much advantage from a roulette here
	// because the size of a species is relatively small.
	// But you can use it by using the roulette code here
	//for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
	//  total_fitness+=(*curorg)->fitness;
	//}

	
	//Check for a mistake
	if ((expected_offspring>0)&&
		(organisms.size()==0)) {
			    std::cout<<"NOT AN ERROR:  ATTEMPT TO REPRODUCE OUT OF EMPTY SPECIES"<<std::endl;
                            exit(0);
			return false;
		}

                                std::cout<<"poolsize: "<<poolsize<<std::endl;

		thechamp=(*(organisms.begin()));

		//Create the designated number of offspring for the Species
		//one at a time
		for (count=0;count<expected_offspring;count++) {

			mut_struct_baby=false;
			mate_baby=false;

			outside=false;

			//Debug Trap
			if (expected_offspring>NEAT::pop_size) {
				//      std::cout<<"ALERT: EXPECTED OFFSPRING = "<<expected_offspring<<std::endl;
				//      cin>>pause;
			}

			//If we have a super_champ (Population champion), finish off some special clones
			//DO NOT NEED ELITIST NSGA-2??
			if ((thechamp->super_champ_offspring) > 0) {
                                std::cout<<"ERROR:: ((thechamp->super_champ_offspring) > 0) not implemented right now"<<std::endl;
                                exit(0);
				mom=thechamp;
				new_genome=(mom->gnome)->duplicate(count);

				if ((thechamp->super_champ_offspring) == 1) {

				}

				//Most superchamp offspring will have their connection weights mutated only
				//The last offspring will be an exact duplicate of this super_champ
				//Note: Superchamp offspring only occur with stolen babies!
				//      Settings used for published experiments did not use this
				if ((thechamp->super_champ_offspring) > 1) {
					if ((randfloat()<0.8)||
						(NEAT::mutate_add_link_prob==0.0)) 
						//ABOVE LINE IS FOR:
						//Make sure no links get added when the system has link adding disabled
						new_genome->mutate_link_weights(mut_power,1.0,GAUSSIAN);
					else {
						//Sometimes we add a link to a superchamp
						net_analogue=new_genome->genesis(generation);
						new_genome->mutate_add_link(pop->innovations,pop->cur_innov_num,NEAT::newlink_tries);
						delete net_analogue;
						mut_struct_baby=true;
					}
				}

				baby=new Organism(0.0, 0.0, new_genome,generation);

				if ((thechamp->super_champ_offspring) == 1) {
					if (thechamp->pop_champ) {
						//std::cout<<"The new org baby's genome is "<<baby->gnome<<std::endl;
						baby->pop_champ_child=true;
						baby->high_fit=mom->orig_fitness1;
					}
				}

				thechamp->super_champ_offspring--;
			}
			//DO NOT NEED TO CLONE THE CHAMP IN THE ELITIST NSGA-2
                        //If we have a Species champion, just clone it  
			else if ((!champ_done)&&
				(expected_offspring>5)) {
                                        std::cout<<"ERROR:: Champ duplication is not required in elitist NSGA-2"<<std::endl;
                                        exit(0);

					mom=thechamp; //Mom is the champ

					new_genome=(mom->gnome)->duplicate(count);

					baby=new Organism(0.0, 0.0, new_genome,generation);  //Baby is just like mommy

					champ_done=true;

				}
				//First, decide whether to mate or mutate
				//If there is only one organism in the pool, then always mutate
			else if ((randfloat()<NEAT::mutate_only_prob)||
				poolsize== 1) {

					//Choose the random parent

					////Roulette Wheel
					//marble=randfloat()*total_fitness;
					//curorg=organisms.begin();
					//spin=(*curorg)->fitness;
					//while(spin<marble) {
					//++curorg;

					////Keep the wheel spinning
					//spin+=(*curorg)->fitness;
					//}
					////Finished roulette
					//
                                        mom = binary_tournament_select(poolsize);

					new_genome=(mom->gnome)->duplicate(count);

					//Do the mutation depending on probabilities of 
					//various mutations

					if (randfloat()<NEAT::mutate_add_node_prob) {
			                        Genome *last_genome=(new_genome)->duplicate(1);
						//std::cout<<"mutate add node"<<std::endl;
						new_genome->mutate_add_node(pop->innovations,pop->cur_node_id,pop->cur_innov_num);
						mut_struct_baby=true;
                                                char *s = "mutate_add_node";
                                                check_lstm_genes(last_genome, new_genome, s);
                                                delete last_genome;
					}
		                 	else if (randfloat()<NEAT::mutate_add_lstm_node_prob) {//Adding LSTM node is less preferred than regular node
			                        Genome *last_genome=(new_genome)->duplicate(1);
		                 		new_genome->mutate_add_lstm_node(pop->innovations,pop->cur_node_id,pop->cur_innov_num);
		                 		mut_struct_baby=true;
                                                char *s = "mutate_add_lstm_node";
                                                check_lstm_genes(last_genome, new_genome, s);
                                                delete last_genome;
                                         }
					else if (randfloat()<NEAT::mutate_add_link_prob) {
						//std::cout<<"mutate add link"<<std::endl;
			                        Genome *last_genome=(new_genome)->duplicate(1);
						net_analogue=new_genome->genesis(generation);
						new_genome->mutate_add_link(pop->innovations,pop->cur_innov_num,NEAT::newlink_tries);
						delete net_analogue;
						mut_struct_baby=true;
                                                char *s = "mutate_add_link";
                                                check_lstm_genes(last_genome, new_genome, s);
                                                delete last_genome;
					}
					//NOTE:  A link CANNOT be added directly after a node was added because the phenotype
					//       will not be appropriately altered to reflect the change
					else {
						//If we didn't do a structural mutation, we do the other kinds

						if (randfloat()<NEAT::mutate_random_trait_prob) {
							//std::cout<<"mutate random trait"<<std::endl;
							new_genome->mutate_random_trait();
						}
						if (randfloat()<NEAT::mutate_link_trait_prob) {
							//std::cout<<"mutate_link_trait"<<std::endl;
							new_genome->mutate_link_trait(1);
						}
						if (randfloat()<NEAT::mutate_node_trait_prob) {
							//std::cout<<"mutate_node_trait"<<std::endl;
							new_genome->mutate_node_trait(1);
						}
						if (randfloat()<NEAT::mutate_link_weights_prob) {
			                                Genome *last_genome=(new_genome)->duplicate(1);
							//std::cout<<"mutate_link_weights"<<std::endl;
							new_genome->mutate_link_weights(mut_power,1.0,GAUSSIAN);
                                                        char *s = "mutate_link_weights";
                                                        check_lstm_genes(last_genome, new_genome, s);
                                                        delete last_genome;
						}
						if (randfloat()<NEAT::mutate_toggle_enable_prob) {
			                                Genome *last_genome=(new_genome)->duplicate(1);
							//std::cout<<"mutate toggle enable"<<std::endl;
							new_genome->mutate_toggle_enable(1);
                                                        char *s = "mutate_toggle_enable";
                                                        check_lstm_genes(last_genome, new_genome, s);
                                                        delete last_genome;

						}
						if (randfloat()<NEAT::mutate_gene_reenable_prob) {
			                                Genome *last_genome=(new_genome)->duplicate(1);
							//std::cout<<"mutate gene reenable"<<std::endl;
							new_genome->mutate_gene_reenable();
                                                        char *s = "mutate_gene_reenable";
                                                        check_lstm_genes(last_genome, new_genome, s);
                                                        delete last_genome;
						}
					}

					baby=new Organism(0.0, 0.0, new_genome,generation);

				}

				//Otherwise we should mate 
			else {

		         	//Choose parents
                                mom = binary_tournament_select(poolsize);

				//Choose random dad

				if ((randfloat()>NEAT::interspecies_mate_rate)) {
					//Mate within Species

                                        dad = binary_tournament_select(poolsize);
				}
				else {

					//Mate outside Species  
					randspecies=this;

					//Select a random species
					giveup=0;  //Give up if you cant find a different Species
					while((randspecies==this)&&(giveup<5)) {

						//This old way just chose any old species
						//randspeciesnum=randint(0,(pop->species).size()-1);

						//Choose a random species tending towards better species
						randmult=gaussrand()/4;
						if (randmult>1.0) randmult=1.0;
						//This tends to select better species
						randspeciesnum=(int) floor((randmult*(sorted_species.size()-1.0))+0.5);
						cursp=(sorted_species.begin());
						for(spcount=0;spcount<randspeciesnum;spcount++)
							++cursp;
						randspecies=(*cursp);

						++giveup;
					}

					//OLD WAY: Choose a random dad from the random species
					//Select a random dad from the random Species
					//NOTE:  It is possible that a mating could take place
					//       here between the mom and a baby from the NEW
					//       generation in some other Species
					//orgnum=randint(0,(randspecies->organisms).size()-1);
					//curorg=(randspecies->organisms).begin();
					//for(orgcount=0;orgcount<orgnum;orgcount++)
					//  ++curorg;
					//dad=(*curorg);            

					//New way: Make dad be a champ from the random species
					dad=(*((randspecies->organisms).begin()));
					

					outside=true;	
				}

				//Perform mating based on probabilities of differrent mating types
				if (randfloat()<NEAT::mate_multipoint_prob) { 
			                Genome *last_genome=(mom->gnome)->duplicate(1);
					new_genome=(mom->gnome)->mate_multipoint(dad->gnome,count,mom->front_num,dad->front_num,mom->crowd_dist,dad->crowd_dist,outside);
                                        char *s = "mate_multipoint";
                                        check_lstm_genes(last_genome, new_genome, s);
                                        delete last_genome;
				}
				else if (randfloat()<(NEAT::mate_multipoint_avg_prob/(NEAT::mate_multipoint_avg_prob+NEAT::mate_singlepoint_prob))) {
			                Genome *last_genome=(mom->gnome)->duplicate(1);
					new_genome=(mom->gnome)->mate_multipoint_avg(dad->gnome,count,mom->front_num,dad->front_num,mom->crowd_dist,dad->crowd_dist,outside);
                                        char *s = "mate_multipoint_avg";
                                        check_lstm_genes(last_genome, new_genome, s);
                                        delete last_genome;
				}
				else {
			                Genome *last_genome=(mom->gnome)->duplicate(1);
					new_genome=(mom->gnome)->mate_singlepoint(dad->gnome,count);
                                        char *s = "mate_singlepoint";
                                        check_lstm_genes(last_genome, new_genome, s);
                                        delete last_genome;
				}

				mate_baby=true;

				//Determine whether to mutate the baby's Genome
				//This is done randomly or if the mom and dad are the same organism
				if ((randfloat()>NEAT::mate_only_prob)||
					((dad->gnome)->genome_id==(mom->gnome)->genome_id)||
					(((dad->gnome)->compatibility(mom->gnome))==0.0))
				{

					//Do the mutation depending on probabilities of 
					//various mutations
					if (randfloat()<NEAT::mutate_add_node_prob) {
			                        Genome *last_genome=(new_genome)->duplicate(1);
						new_genome->mutate_add_node(pop->innovations,pop->cur_node_id,pop->cur_innov_num);
						//  std::cout<<"mutate_add_node: "<<new_genome<<std::endl;
						mut_struct_baby=true;
                                                char *s = "mutate_add_node";
                                                check_lstm_genes(last_genome, new_genome, s);
                                                delete last_genome;
					}
		                 	else if (randfloat()<NEAT::mutate_add_lstm_node_prob) {//Adding LSTM node is less preferred than regular node
			                        Genome *last_genome=(new_genome)->duplicate(1);
		                 		new_genome->mutate_add_lstm_node(pop->innovations,pop->cur_node_id,pop->cur_innov_num);
		                 		mut_struct_baby=true;
                                                char *s = "mutate_add_lstm_node";
                                                check_lstm_genes(last_genome, new_genome, s);
                                                delete last_genome;
                                        }
					else if (randfloat()<NEAT::mutate_add_link_prob) {
			                        Genome *last_genome=(new_genome)->duplicate(1);
						net_analogue=new_genome->genesis(generation);
						new_genome->mutate_add_link(pop->innovations,pop->cur_innov_num,NEAT::newlink_tries);
						delete net_analogue;
						//std::cout<<"mutate_add_link: "<<new_genome<<std::endl;
						mut_struct_baby=true;
                                                char *s = "mutate_add_link";
                                                check_lstm_genes(last_genome, new_genome, s);
                                                delete last_genome;
					}
					else {
						//Only do other mutations when not doing sturctural mutations

						if (randfloat()<NEAT::mutate_random_trait_prob) {
							new_genome->mutate_random_trait();
							//std::cout<<"..mutate random trait: "<<new_genome<<std::endl;
						}
						if (randfloat()<NEAT::mutate_link_trait_prob) {
							new_genome->mutate_link_trait(1);
							//std::cout<<"..mutate link trait: "<<new_genome<<std::endl;
						}
						if (randfloat()<NEAT::mutate_node_trait_prob) {
							new_genome->mutate_node_trait(1);
							//std::cout<<"mutate_node_trait: "<<new_genome<<std::endl;
						}
						if (randfloat()<NEAT::mutate_link_weights_prob) {
			                                Genome *last_genome=(new_genome)->duplicate(1);
							new_genome->mutate_link_weights(mut_power,1.0,GAUSSIAN);
							//std::cout<<"mutate_link_weights: "<<new_genome<<std::endl;
                                                        char *s = "mutate_line_weights";
                                                        check_lstm_genes(last_genome, new_genome, s);
                                                        delete last_genome;
						}
						if (randfloat()<NEAT::mutate_toggle_enable_prob) {
			                                Genome *last_genome=(new_genome)->duplicate(1);
							new_genome->mutate_toggle_enable(1);
							//std::cout<<"mutate_toggle_enable: "<<new_genome<<std::endl;
                                                        char *s = "mutate_toggle_enable";
                                                        check_lstm_genes(last_genome, new_genome, s);
                                                        delete last_genome;
						}
						if (randfloat()<NEAT::mutate_gene_reenable_prob) {
			                                Genome *last_genome=(new_genome)->duplicate(1);
							new_genome->mutate_gene_reenable(); 
							//std::cout<<"mutate_gene_reenable: "<<new_genome<<std::endl;
                                                        char *s = "mutate_gene_reenable";
                                                        check_lstm_genes(last_genome, new_genome, s);
                                                        delete last_genome;
						}
					}

					//Create the baby
					baby=new Organism(0.0, 0.0, new_genome,generation);

				}
				else {
					//Create the baby without mutating first
					baby=new Organism(0.0, 0.0, new_genome,generation);
				}

			}

			//Add the baby to its proper Species
			//If it doesn't fit a Species, create a new one

			baby->mut_struct_baby=mut_struct_baby;
			baby->mate_baby=mate_baby;

			curspecies=(pop->species).begin();
			if (curspecies==(pop->species).end()){
				//Create the first species
				newspecies=new Species(++(pop->last_species),true);
				(pop->species).push_back(newspecies);
				newspecies->add_Organism(baby);  //Add the baby
				baby->species=newspecies;  //Point the baby to its species
			} 
			else {
				comporg=(*curspecies)->first();
				found=false;
				while((curspecies!=(pop->species).end())&&
					(!found)) {	
						if (comporg==0) {
							//Keep searching for a matching species
							++curspecies;
							if (curspecies!=(pop->species).end())
								comporg=(*curspecies)->first();
						}
						else if (((baby->gnome)->compatibility(comporg->gnome))<NEAT::compat_threshold) {
							//Found compatible species, so add this organism to it
							(*curspecies)->add_Organism(baby);
							baby->species=(*curspecies);  //Point organism to its species
							found=true;  //Note the search is over
						}
						else {
							//Keep searching for a matching species
							++curspecies;
							if (curspecies!=(pop->species).end()) 
								comporg=(*curspecies)->first();
						}
					}

					//If we didn't find a match, create a new species
					if (found==false) {
					  newspecies=new Species(++(pop->last_species),true);
					  //std::std::cout<<"CREATING NEW SPECIES "<<pop->last_species<<std::std::endl;
					  (pop->species).push_back(newspecies);
					  newspecies->add_Organism(baby);  //Add the baby
					  baby->species=newspecies;  //Point baby to its species
					}


			} //end else 

		}



		return true;
}

bool NEAT::order_species(Species *x, Species *y) { 
	//std::cout<<"Comparing "<<((*((x->organisms).begin()))->orig_fitness)<<" and "<<((*((y->organisms).begin()))->orig_fitness)<<": "<<(((*((x->organisms).begin()))->orig_fitness) > ((*((y->organisms).begin()))->orig_fitness))<<std::endl;
	return (((*((x->organisms).begin()))->orig_fitness1) > ((*((y->organisms).begin()))->orig_fitness1));
}

bool NEAT::order_species_by_front_num_crowd_dist(Species *x, Species *y) { 
	//std::cout<<"Comparing "<<((*((x->organisms).begin()))->orig_fitness)<<" and "<<((*((y->organisms).begin()))->orig_fitness)<<": "<<(((*((x->organisms).begin()))->orig_fitness) > ((*((y->organisms).begin()))->orig_fitness))<<std::endl;
	return ((((*((x->organisms).begin()))->front_num) <  ((*((y->organisms).begin()))->front_num)) ||
               ((((*((x->organisms).begin()))->front_num) == ((*((y->organisms).begin()))->front_num)) && 
                    (((*((x->organisms).begin()))->crowd_dist) > ((*((y->organisms).begin()))->crowd_dist))));
}

bool NEAT::order_new_species(Species *x, Species *y) {
	return (x->compute_max_fitness() > y->compute_max_fitness());
}


