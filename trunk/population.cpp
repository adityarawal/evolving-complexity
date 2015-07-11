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
#include "population.h"
#include "organism.h"
#include <iostream>
#include <sstream>
#include <fstream>
using namespace NEAT;

Population::Population(Genome *g,int size) {
	winnergen=0;
	highest_fitness1=0.0;
	highest_last_changed1=0;
	highest_fitness2=0.0;
	highest_last_changed2=0;
	spawn(g,size);
}

Population::Population(Genome *g,int size, float power) {
	winnergen=0;
	highest_fitness1=0.0;
	highest_last_changed1=0;
	highest_fitness2=0.0;
	highest_last_changed2=0;
	clone(g, size, power);
}

//Population::Population(int size,int i,int o, int nmax, bool r, double linkprob) {    
//int count;
//Genome *new_genome; 

//cout<<"Making a random pop"<<endl;

//winnergen=0;
//highest_fitness=0.0;
//highest_last_changed=0;

//for(count=0;count<size;count++) {
//new_genome=new Genome(count,i,o,randint(0,nmax),nmax,r,linkprob);
//organisms.push_back(new Organism(0,new_genome,1));
//}

//cur_node_id=i+o+nmax+1;;
//cur_innov_num=(i+o+nmax)*(i+o+nmax)+1;

//cout<<"Calling speciate"<<endl;
//speciate(); 

//}


//MSC Addition
//Added the ability for a population to be spawned
//off of a vector of Genomes.  Useful when converging.
Population::Population(std::vector<Genome*> genomeList, float power) {
	
	winnergen=0;
	highest_fitness1=0.0;
	highest_last_changed1=0;
	highest_fitness2=0.0;
	highest_last_changed2=0;
		
	int count;
	Genome *new_genome;
	Organism *new_organism;

	//Create size copies of the Genome
	//Start with perturbed linkweights
	for (std::vector<Genome*>::iterator iter = genomeList.begin(); iter != genomeList.end(); ++iter)
	{

		new_genome=(*iter); 
		if(power>0)
			new_genome->mutate_link_weights(power,1.0,GAUSSIAN);
		//new_genome->mutate_link_weights(1.0,1.0,COLDGAUSSIAN);
		new_genome->randomize_traits();
		new_organism=new Organism(0.0,0.0, new_genome,1);
		organisms.push_back(new_organism);
	}

	//Keep a record of the innovation and node number we are on
	cur_node_id=new_genome->get_last_node_id();
	cur_innov_num=new_genome->get_last_gene_innovnum();

	//Separate the new Population into species
	speciate();
}

Population::Population(const char *filename) {

	char curword[128];  //max word size of 128 characters
	char curline[1024]; //max line size of 1024 characters
	char delimiters[] = " \n";

	Genome *new_genome;

	winnergen=0;

	highest_fitness1=0.0;
	highest_last_changed1=0;
	highest_fitness2=0.0;
	highest_last_changed2=0;

	cur_node_id=0;
	cur_innov_num=0.0;

	int curwordnum = 0;

	std::ifstream iFile(filename);
	if (!iFile) {
		printf("Can't open genomes file for input");
		return;
	}

	else {
		bool md = false;
		char metadata[128];
		//Loop until file is finished, parsing each line
		while (!iFile.eof()) 
		{
			iFile.getline(curline, sizeof(curline));
            std::stringstream ss(curline);
			//strcpy(curword, NEAT::getUnit(curline, 0, delimiters));
            ss >> curword;
            //std::cout << curline << std::endl;

			//Check for next
			if (strcmp(curword,"genomestart")==0) 
			{
				//strcpy(curword, NEAT::getUnit(curline, 1, delimiters));
				//int idcheck = atoi(curword);

                int idcheck;
                ss >> idcheck;

				// If there isn't metadata, set metadata to ""
				if(md == false)  {
					strcpy(metadata, "");
				}
				md = false;

				new_genome=new Genome(idcheck,iFile);
				organisms.push_back(new Organism(0.0, 0.0, new_genome,1, metadata));
				if (cur_node_id<(new_genome->get_last_node_id()))
					cur_node_id=new_genome->get_last_node_id();

				if (cur_innov_num<(new_genome->get_last_gene_innovnum()))
					cur_innov_num=new_genome->get_last_gene_innovnum();
			}
			else if (strcmp(curword,"/*")==0) 
			{
				// New metadata possibly, so clear out the metadata
				strcpy(metadata, "");
				curwordnum=1;
				//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
                ss >> curword;

				while(strcmp(curword,"*/")!=0)
				{
					// If we've started to form the metadata, put a space in the front
					if(md) {
						strncat(metadata, " ", 128 - strlen(metadata));
					}

					// Append the next word to the metadata, and say that there is metadata
					strncat(metadata, curword, 128 - strlen(metadata));
					md = true;

					//strcpy(curword, NEAT::getUnit(curline, curwordnum++, delimiters));
                    ss >> curword;
				}
			}
			//Ignore comments - they get printed to screen
			//else if (strcmp(curword,"/*")==0) {
			//	iFile>>curword;
			//	while (strcmp(curword,"*/")!=0) {
			//cout<<curword<<" ";
			//		iFile>>curword;
			//	}
			//	cout<<endl;

			//}
			//Ignore comments surrounded by - they get printed to screen
		}

		iFile.close();

		speciate();

	}
}


Population::~Population() {

	std::vector<Species*>::iterator curspec;
	std::vector<Organism*>::iterator curorg;
	//std::vector<Generation_viz*>::iterator cursnap;

	if (species.begin()!=species.end()) {
		for(curspec=species.begin();curspec!=species.end();++curspec) {
			delete (*curspec);
		}
	}
	else {
		for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
			delete (*curorg);
		}
	}

	for (std::vector<Innovation*>::iterator iter = innovations.begin(); iter != innovations.end(); ++iter)
		delete *iter;

	//Delete the snapshots
	//		for(cursnap=generation_snapshots.begin();cursnap!=generation_snapshots.end();++cursnap) {
	//			delete (*cursnap);
	//		}
}

bool Population::verify() {
	std::vector<Organism*>::iterator curorg;

	bool verification;

	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		verification=((*curorg)->gnome)->verify();
	}

	return verification;
} 

bool Population::clone(Genome *g,int size, float power) {
	int count;
	Genome *new_genome;
	Organism *new_organism;

	new_genome = g->duplicate(1); 
	new_organism = new Organism(0.0, 0.0, new_genome,1);
	organisms.push_back(new_organism);
	
	//Create size copies of the Genome
	//Start with perturbed linkweights
	for(count=2;count<=size;count++) {
		//cout<<"CREATING ORGANISM "<<count<<endl;
		new_genome=g->duplicate(count); 
		if(power>0)
			new_genome->mutate_link_weights(power,1.0,GAUSSIAN);
		
		new_genome->randomize_traits();
		new_organism=new Organism(0.0, 0.0, new_genome,1);
		organisms.push_back(new_organism);
	}

	//Keep a record of the innovation and node number we are on
	cur_node_id=new_genome->get_last_node_id();
	cur_innov_num=new_genome->get_last_gene_innovnum();

	//Separate the new Population into species
	speciate();

	return true;
}

bool Population::spawn(Genome *g,int size) {
	int count;
	Genome *new_genome;
	Organism *new_organism;

	//Create size copies of the Genome
	//Start with perturbed linkweights
	for(count=1;count<=size;count++) {
		//cout<<"CREATING ORGANISM "<<count<<endl;

		new_genome=g->duplicate(count); 
		//new_genome->mutate_link_weights(1.0,1.0,GAUSSIAN);
		//new_genome->mutate_link_weights(1.0,1.0,COLDGAUSSIAN);
		//new_genome->randomize_traits(); Aditya: No need to change traits since they are not being used
		new_organism=new Organism(0.0, 0.0, new_genome,1);//(Aditya: for NSGA-2)
		organisms.push_back(new_organism);
	}

	//Keep a record of the innovation and node number we are on
	cur_node_id=new_genome->get_last_node_id();
	cur_innov_num=new_genome->get_last_gene_innovnum();

	//Separate the new Population into species
	speciate();

	return true;
}

bool Population::speciate() {
	std::vector<Organism*>::iterator curorg;  //For stepping through Population
	std::vector<Species*>::iterator curspecies; //Steps through species
	Organism *comporg=0;  //Organism for comparison 
	Species *newspecies; //For adding a new species

	int counter=0; //Species counter

	//Step through all existing organisms
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {

		//For each organism, search for a species it is compatible to
		curspecies=species.begin();
		if (curspecies==species.end()){
			//Create the first species
			newspecies=new Species(++counter);
			species.push_back(newspecies);
			newspecies->add_Organism(*curorg);  //Add the current organism
			(*curorg)->species=newspecies;  //Point organism to its species
		} 
		else {
			comporg=(*curspecies)->first();
			while((comporg!=0)&&
				(curspecies!=species.end())) {

					if ((((*curorg)->gnome)->compatibility(comporg->gnome))<NEAT::compat_threshold) {

						//Found compatible species, so add this organism to it
						(*curspecies)->add_Organism(*curorg);
						(*curorg)->species=(*curspecies);  //Point organism to its species
						comporg=0;  //Note the search is over
					}
					else {

						//Keep searching for a matching species
						++curspecies;
						if (curspecies!=species.end()) 
							comporg=(*curspecies)->first();
					}
				}

				//If we didn't find a match, create a new species
				if (comporg!=0) {
					newspecies=new Species(++counter);
					species.push_back(newspecies);
					newspecies->add_Organism(*curorg);  //Add the current organism
					(*curorg)->species=newspecies;  //Point organism to its species
				}

		} //end else 

	} //end for

	last_species=counter;  //Keep track of highest species

	return true;
}

bool Population::print_to_file_by_species(char *filename) {

  std::vector<Species*>::iterator curspecies;

  std::ofstream outFile(filename,std::ios::out);

  //Make sure it worked
  if (!outFile) {
    std::cerr<<"Can't open "<<filename<<" for output"<<std::endl;
    return false;
  }


  //Step through the Species and print them to the file
  for(curspecies=species.begin();curspecies!=species.end();++curspecies) {
    (*curspecies)->print_to_file(outFile);
  }

  outFile.close();

  return true;

}


bool Population::print_to_file_by_species(std::ostream& outFile) {

	std::vector<Species*>::iterator curspecies;

	//ofstream outFile(filename,ios::out);
	//std::ostream outFile;
	//ResourceManager->openFileForWrite(outFile, fileName, std::ostream::Write);

	//Make sure it worked
	//if (!outFile) {
	//	cerr<<"Can't open "<<filename<<" for output"<<endl;
	//	return false;
	//}


	//Step through the Species and print them to the file
	for(curspecies=species.begin();curspecies!=species.end();++curspecies) {
		(*curspecies)->print_to_file(outFile);
	}

	return true;

}

bool NEAT::order_org_fitness1(Organism *x, Organism *y) { 
	return (x->fitness1 > y->fitness1); //Sort in descending order of fitness 1
}

bool NEAT::order_org_fitness2(Organism *x, Organism *y) { 
	return (x->fitness2 > y->fitness2); //Sort in descending order of fitness 2
}

bool NEAT::order_org_crowd_dist(Organism *x, Organism *y) { 
	return (x->crowd_dist > y->crowd_dist); //Sort in descending order of crowding distance
}

bool Population::p_dominates_q (Organism* p, Organism* q) {
        if ((p->fitness1 < q->fitness1) || (p->fitness2 < q->fitness2)) {
                return false;
        }
        else if ((p->fitness1 > q->fitness1) || (p->fitness2 > q->fitness2)) {
                return true;
        }
        else if ((p->fitness1 == q->fitness1) && (p->fitness2 == q->fitness2)) {
                return false;
        }
        else {
                std::cout <<"ERROR in p_dominates_q()... Exiting"<<std::endl;
                exit(0);
        }
}

std::vector <double> Population::compute_fitness_range() {

	std::vector<Organism*>::iterator curorg;
	std::vector<Organism*> sorted_organisms;  //Organisms sorted by fitness
        std::vector <double> fitness_range;
	
        //Stick the Organism pointers into a new organism list for sorting
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		sorted_organisms.push_back(*curorg);
	}

        double temp_range; 
        //Sort the Organisms by fitness objective 1
        std::sort(sorted_organisms.begin(), sorted_organisms.end(), order_org_fitness1);
        //Compute absolute difference in max and min fitness values
        temp_range = std::abs((*(sorted_organisms.begin()))->fitness1 - (*(--sorted_organisms.end()))->fitness1);
        fitness_range.push_back(temp_range);

        //Sort the Organisms by fitness objective 2
        std::sort(sorted_organisms.begin(), sorted_organisms.end(), order_org_fitness2);
        //Compute absolute difference in max and min fitness values
        temp_range = std::abs((*(sorted_organisms.begin()))->fitness2 - (*(--sorted_organisms.end()))->fitness2);
        fitness_range.push_back(temp_range);
        return fitness_range; 
}

std::vector< std::vector<Organism*> > Population::assign_domination_fronts(){
        std::vector<Organism*>::iterator p;
        std::vector<Organism*>::iterator q;
	std::vector<Organism*> curfront;
	std::vector< std::vector<Organism*> > fronts;

        for(p=organisms.begin();p!=organisms.end();++p) {
                (*p)->domination_list.clear(); //Reset parameters
                (*p)->num_dominant = 0;        //Reset parameters
                (*p)->front_num = -1;//Most dominant front
                for(q=organisms.begin();q!=organisms.end();++q) {
                        if (p_dominates_q((*p),(*q))){//If p dominates q (Dereferencing iterator to get pointer)
                                (*p)->domination_list.push_back((*q));//S_p = S_p U {q}
                        }
                        else if (p_dominates_q((*q),(*p))) { //If q dominates p
                                (*p)->num_dominant = (*p)->num_dominant + 1;//n_p = n_p + 1
                        }
                }
                if ((*p)->num_dominant == 0) {
                        (*p)->front_num = 0;//Most dominant front
                        curfront.push_back((*p));
                }
	}
        fronts.push_back(curfront); //Added the first/most-dominant front
        
        int i = 0;
        while (true) {
                curfront.clear();//Q = {}
                for(p=fronts[i].begin();p!=fronts[i].end();++p) {//for each p in front_i
                        for(q=(*p)->domination_list.begin();q!=(*p)->domination_list.end();++q) {//for each q in S_p
                                (*q)->num_dominant = (*q)->num_dominant-1;//n_q = n_q - 1
                                if ((*q)->num_dominant == 0) {//If this is the next dominant front
                                        (*q)->front_num = i+1; //q_rank = i+1
                                        curfront.push_back((*q)); //Q = Q U {q}
                                }
                        }
                }
                i=i+1;
                if (!curfront.empty()) {//Add front only if it not empty
                        fronts.push_back(curfront);//front_i = Q
                }
                else {
                        break;//No more fronts remaining
                } 
        }
        return fronts;
}

//Assign crowding distance within front
void Population::assign_crowding_distance(std::vector<Organism*> front, int num_obj) {
        double fitness_range;
        int length = front.size();
        std::vector<Organism*>::iterator curorg;

        //Initialize crowding distance to zero
        for(curorg=front.begin();curorg!=front.end();++curorg) {
                (*curorg)->crowd_dist = 0.0;
        }

        int obj_count = 0;
        while (obj_count < num_obj) {
               //Compute Crowding distance for fitness objectives
               if (obj_count == 0) {
                       //Sort the front by fitness objective 1
                       std::sort(front.begin(), front.end(), order_org_fitness1);
                       fitness_range =  (*(front.begin()))->fitness1 - (*(--front.end()))->fitness1;
               }
               else if (obj_count == 1) {
                       //Sort the front by fitness objective 1
                       std::sort(front.begin(), front.end(), order_org_fitness2);
                       fitness_range =  (*(front.begin()))->fitness2 - (*(--front.end()))->fitness2;
               }

               //First and last organisms are always selected  
               (*(front.begin()))->crowd_dist = 100000.0; //Infinity 
               (*(--front.end()))->crowd_dist = 100000.0; //Infinity

               double temp_dist; 
               for (int i=1; i<length-1; i++) {//Start from the second organism and end at second last organism
                       if (fitness_range <= 0) {
                               std::cout << "ERROR:: Fitness Range cannot be less than equal to ZERO"<<std::endl;
                               fitness_range = 1;//Is this CORRECT?? 
                               //exit(0);
                       }
                       if (obj_count==0) {//For objective 1
                               temp_dist = std::abs(front[i+1]->fitness1 - front[i-1]->fitness1);
                       }
                       else if (obj_count==1){//For objective 2
                               temp_dist = std::abs(front[i+1]->fitness2 - front[i-1]->fitness2);
                       }
                       front[i]->crowd_dist = front[i]->crowd_dist + (temp_dist/fitness_range); 
               }
               obj_count = obj_count + 1;
        }
}

bool Population::epoch_multiobj(int generation) {

        std::vector<Species*>::iterator curspecies;

	std::vector<Organism*>::iterator curorg;
	std::vector<Organism*>::iterator deadorg;

	std::vector<Innovation*>::iterator curinnov;  
	std::vector<Innovation*>::iterator deadinnov;  //For removing old Innovs

	double total=0.0; //Used to compute average fitness over all Organisms

	double overall_average1;  //The average modified fitness among ALL organisms
	double overall_average2;  //The average modified fitness among ALL organisms

	int orgcount;

	//The fractional parts of expected offspring that can be 
	//Used only when they accumulate above 1 for the purposes of counting
	//Offspring
	double skim; 
	int total_expected;  //precision checking
	int total_organisms=organisms.size();
	int max_expected;
	Species *best_species;
	int final_expected;

	int pause;

	int half_pop;

	int best_species_num;  //Used in debugging to see why (if) best species dies
	bool best_ok;

	//We can try to keep the number of species constant at this number
	int num_species=species.size();

	std::vector< std::vector<Organism*> > nondominated_fronts;
	std::vector<Organism*> cur_front;
        nondominated_fronts = assign_domination_fronts();
        int num_obj = 2; //Number of objectives in the problem
       
        //Compute fitness range for each objective
        //fitness_range = compute_fitness_range();

	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
                (*curorg)->eliminate = true; //Mark all the organisms for elimination initiallly
	}

        //Go through all the fronts and find the top organisms (non-dominated) pop_size individuals
        int count_org_selected = 0;
        int curfront = 0;
        while ((nondominated_fronts[curfront].size() + count_org_selected) <= (NEAT::pop_size/2)) {//Only top N are selected for parenting (pop_size = 2N)
	        for(curorg=nondominated_fronts[curfront].begin();curorg!=nondominated_fronts[curfront].end();++curorg) {
                        (*curorg)->eliminate = false; //Selected organisms are not eliminated 
	        }
                assign_crowding_distance(nondominated_fronts[curfront], num_obj); //Assign crowding distance within front
                count_org_selected = count_org_selected + nondominated_fronts[curfront].size();
                curfront = curfront+1;
        }
        //If selecting a front causes population overflow
        //Select subset of organisms from this front 
        if (count_org_selected < (NEAT::pop_size)/2) {
                if (nondominated_fronts[curfront].empty()) {//For DEBUG purpose, add this check
                        std::cout<<"ERRORR in NSGA-2 front assignment"<<std::endl;
                        exit(0);
                } 
                int count_org_remaining;
                count_org_remaining = (NEAT::pop_size/2)-count_org_selected;
                assign_crowding_distance(nondominated_fronts[curfront], num_obj);

                //Sort nondominated_fronts[curfront] according to the descending order of crowding distance 
                std::sort(nondominated_fronts[curfront].begin(), nondominated_fronts[curfront].end(), order_org_crowd_dist);
                //Pick top count_org_remaining organisms from the sorted front (set eliminate to false)
	        for(curorg=nondominated_fronts[curfront].begin();curorg!=nondominated_fronts[curfront].begin() + count_org_remaining;++curorg) {
                        (*curorg)->eliminate = false; //Selected organisms are not eliminated (CHECK THIS) 
	        }
        }
	//Kill off all Organisms marked for death.  The remainder
	//will be allowed to reproduce.
	curorg=organisms.begin();
	while(curorg!=organisms.end()) {
		if (((*curorg)->eliminate)) {
			//Remove the organism from its Species
			((*curorg)->species)->remove_org(*curorg);

			//Delete the organism from memory
			delete (*curorg);

			//Remember where we are
			deadorg=curorg;
			++curorg;

			//iter2 =  v.erase(iter); 

			//Remove the organism from the master list
			curorg=organisms.erase(deadorg);

		}
		else {
			++curorg;
		}

	}

	std::cout<<"Number of Species: "<<num_species<<std::endl;

	//Go through the organisms and add up their fitnesses to compute the
	//overall average
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		total+=(*curorg)->fitness1;
	}
	overall_average1=total/total_organisms;
        total = 0.0;
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		total+=(*curorg)->fitness2;
	}
	overall_average2=total/total_organisms;
        //std::ofstream file1 ("average_dump.txt", std::ios::app);
        //file1 << overall_average1<<" "<<overall_average2<< std::endl;
        //file1.close();
        
	//Check for Population-level stagnation
        double cur_highest_fitness1 = 0.0;
        double cur_highest_fitness2 = 0.0;
        int index = 0; 
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
                if ((*curorg)->fitness1 > cur_highest_fitness1) {
                        cur_highest_fitness1 = (*curorg)->fitness1;
                }
                index = index + 1;
        }
        organisms[index-1]->pop_champ=true;//DEBUG marker of the best of pop (VERIFY THIS)
	index = 0;
        for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
                if ((*curorg)->fitness2 > cur_highest_fitness2) {
                        cur_highest_fitness2 = (*curorg)->fitness2;
                }
                index = index + 1;
        }
	std::cout<<"Generation "<<generation<<": "<<"overall_average 1 = "<<overall_average1<<"  overall_average 2 = "<<overall_average2<<"  Current_highest_fitness1 = "<<cur_highest_fitness1<<"  Current_highest_fitness2 = "<<cur_highest_fitness2<<std::endl;
        
	if (cur_highest_fitness1 > highest_fitness1) {
			highest_fitness1=cur_highest_fitness1;
			highest_last_changed1=0;
			std::cout<<"NEW POPULATION RECORD FITNESS 1: "<<highest_fitness1<<std::endl;
		}
	else {
		++highest_last_changed1;
		std::cout<<highest_last_changed1<<" generations since last population fitness1 record: "<<highest_fitness1<<std::endl;
	}
	if (cur_highest_fitness2 > highest_fitness2) {
			highest_fitness2=cur_highest_fitness2;
			highest_last_changed2=0;
			std::cout<<"NEW POPULATION RECORD FITNESS 2: "<<highest_fitness2<<std::endl;
		}
	else {
		++highest_last_changed2;
		std::cout<<highest_last_changed2<<" generations since last population fitness2 record: "<<highest_fitness2<<std::endl;
	}



        //Reproduce
	curspecies=species.begin();//Only one species in NSGA-2
        (*curspecies)->reproduce_multiobj(generation, this);
        //(*curspecies)->reproduce(generation, this, sorted_species);
        
        //Clearing organisms to keep the code behavior similar to original NEAT code
        //After clearing, the master org list is re-populated from the species list
        //with new genome_id starting with zero.
        //Not Deleting the organism itself because NSGA-2 is elitist and we need to 
        //keep these N (pop_size/2) organisms
        //NOT SURE IF THIS IS CORRECT
        organisms.clear(); 
	
	//Go through the organisms of the curspecies and add them to 
	//the master list
	orgcount=0;
        //++((*curspecies)->age);
	//(*curspecies)->novel=false;
	for(curorg=((*curspecies)->organisms).begin();curorg!=((*curspecies)->organisms).end();++curorg) {
		((*curorg)->gnome)->genome_id=orgcount++;//Reset the genome_id (Mostly used for display purpose
		organisms.push_back(*curorg);
	}

        //Verify that the org in species is the same as that in population
        if ( ((*curspecies)->organisms).size() != organisms.size()) {
                std::cout<<"ERROR: Number of organism mismatch between species and population"<<std::endl;
                exit(0);
        }
	
	//Remove the innovations of the current generation
	curinnov=innovations.begin();
	while(curinnov!=innovations.end()) {
		delete (*curinnov);

		deadinnov=curinnov;
		++curinnov;

		curinnov=innovations.erase(deadinnov);
	}

	//cout<<"Epoch complete"<<endl; 

	return true;

}

bool Population::epoch(int generation) {

	std::vector<Species*>::iterator curspecies;
	std::vector<Species*>::iterator deadspecies;  //For removing empty Species

	std::vector<Organism*>::iterator curorg;
	std::vector<Organism*>::iterator deadorg;

	std::vector<Innovation*>::iterator curinnov;  
	std::vector<Innovation*>::iterator deadinnov;  //For removing old Innovs

	double total=0.0; //Used to compute average fitness over all Organisms

	double overall_average1;  //The average modified fitness among ALL organisms
	double overall_average2;  //The average modified fitness among ALL organisms

	int orgcount;

	//The fractional parts of expected offspring that can be 
	//Used only when they accumulate above 1 for the purposes of counting
	//Offspring
	double skim; 
	int total_expected;  //precision checking
	int total_organisms=organisms.size();
	int max_expected;
	Species *best_species;
	int final_expected;

	int pause;

	//Rights to make babies can be stolen from inferior species
	//and given to their superiors, in order to concentrate exploration on
	//the best species
	int NUM_STOLEN=NEAT::babies_stolen; //Number of babies to steal
	int one_fifth_stolen;
	int one_tenth_stolen;

	std::vector<Species*> sorted_species;  //Species sorted by max fit org in Species
	int stolen_babies; //Babies taken from the bad species and given to the champs

	int half_pop;

	int best_species_num;  //Used in debugging to see why (if) best species dies
	bool best_ok;

	//We can try to keep the number of species constant at this number
	int num_species_target=4;
	int num_species=species.size();
	double compat_mod=0.3;  //Modify compat thresh to control speciation

	std::vector< std::vector<Organism*> > nondominated_fronts;
	std::vector<Organism*> cur_front;
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		cur_front.push_back(*curorg);
	}
        nondominated_fronts.push_back(cur_front);
        nondominated_fronts.push_back(cur_front);
        std::cout<<"Aditya: DONE"<<std::endl;
        exit(0);
        

	//Keeping species diverse
	//This commented out code forces the system to aim for 
	// num_species species at all times, enforcing diversity
	//This tinkers with the compatibility threshold, which
	// normally would be held constant
	/*
	if (generation>1) {
		if (num_species<num_species_target)
			NEAT::compat_threshold-=compat_mod;
		else if (num_species>num_species_target)
			NEAT::compat_threshold+=compat_mod;

		if (NEAT::compat_threshold<0.3) NEAT::compat_threshold=0.3;

	}
	*/


	//Stick the Species pointers into a new Species list for sorting
	for(curspecies=species.begin();curspecies!=species.end();++curspecies) {
		sorted_species.push_back(*curspecies);
	}

	//Sort the Species by max fitness (Use an extra list to do this)
	//These need to use ORIGINAL fitness
	//sorted_species.qsort(order_species);
    std::sort(sorted_species.begin(), sorted_species.end(), order_species);

	//Flag the lowest performing species over age 20 every 30 generations 
	//NOTE: THIS IS FOR COMPETITIVE COEVOLUTION STAGNATION DETECTION

	curspecies=sorted_species.end();
	curspecies--;
	while((curspecies!=sorted_species.begin())&&
		((*curspecies)->age<20))
		--curspecies;
	if ((generation%30)==0)
        {     
                (*curspecies)->obliterate=true;
        }


	std::cout<<"Number of Species: "<<num_species<<std::endl;
	std::cout<<"compat_thresh: "<<compat_threshold<<std::endl;

	//Use Species' ages to modify the objective fitness of organisms
	// in other words, make it more fair for younger species
	// so they have a chance to take hold
	//Also penalize stagnant species
	//Then adjust the fitness using the species size to "share" fitness
	//within a species.
	//Then, within each Species, mark for death 
	//those below survival_thresh*average
	for(curspecies=species.begin();curspecies!=species.end();++curspecies) {
		(*curspecies)->adjust_fitness(); //Aditya: Do not require with NSGA-2
	}

	//Go through the organisms and add up their fitnesses to compute the
	//overall average
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		total+=(*curorg)->fitness1;
	}
	overall_average1=total/total_organisms;
        total = 0.0;
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		total+=(*curorg)->fitness2;
	}
	overall_average2=total/total_organisms;
	std::cout<<"Generation "<<generation<<": "<<"overall_average 1 = "<<overall_average1<<"overall_average 2 = "<<overall_average2<<std::endl;

        //*******START COMMENT (Aditya - NSGA 2)************************* 
	////Now compute expected number of offspring for each individual organism
	//for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
	//	(*curorg)->expected_offspring=(((*curorg)->fitness1)/overall_average);
	//}

	////Now add those offspring up within each Species to get the number of
	////offspring per Species
	//skim=0.0;
	//total_expected=0;
	//for(curspecies=species.begin();curspecies!=species.end();++curspecies) {
	//	skim=(*curspecies)->count_offspring(skim);
	//	total_expected+=(*curspecies)->expected_offspring;
	//}    

	////Need to make up for lost foating point precision in offspring assignment
	////If we lost precision, give an extra baby to the best Species
	//if (total_expected<total_organisms) {
	//	//Find the Species expecting the most
	//	max_expected=0;
	//	final_expected=0;
	//	for(curspecies=species.begin();curspecies!=species.end();++curspecies) {
	//		if ((*curspecies)->expected_offspring>=max_expected) {
	//			max_expected=(*curspecies)->expected_offspring;
	//			best_species=(*curspecies);
	//		}
	//		final_expected+=(*curspecies)->expected_offspring;
	//	}
	//	//Give the extra offspring to the best species
	//	++(best_species->expected_offspring);
	//	final_expected++;

	//	//If we still arent at total, there is a problem
	//	//Note that this can happen if a stagnant Species
	//	//dominates the population and then gets killed off by its age
	//	//Then the whole population plummets in fitness
	//	//If the average fitness is allowed to hit 0, then we no longer have 
	//	//an average we can use to assign offspring.
	//	if (final_expected<total_organisms) {
	//		//      cout<<"Population died!"<<endl;
	//		//cin>>pause;
	//		for(curspecies=species.begin();curspecies!=species.end();++curspecies) {
	//			(*curspecies)->expected_offspring=0;
	//		}
	//		best_species->expected_offspring=total_organisms;
	//	}
	//}
        //*******END COMMENT (Aditya - NSGA 2)************************* 

	//Sort the Species by max fitness (Use an extra list to do this)
	//These need to use ORIGINAL fitness
	//sorted_species.qsort(order_species);
    std::sort(sorted_species.begin(), sorted_species.end(), order_species);

	best_species_num=(*(sorted_species.begin()))->id;

	for(curspecies=sorted_species.begin();curspecies!=sorted_species.end();++curspecies) {

		//Print out for Debugging/viewing what's going on 
		std::cout<<"orig fitness of Species"<<(*curspecies)->id<<"(Size "<<(*curspecies)->organisms.size()<<"): "<<(*((*curspecies)->organisms).begin())->orig_fitness<<" last improved "<<((*curspecies)->age-(*curspecies)->age_of_last_improvement)<<std::endl;
	}

	//Check for Population-level stagnation
	curspecies=sorted_species.begin();
	(*(((*curspecies)->organisms).begin()))->pop_champ=true; //DEBUG marker of the best of pop
	if (((*(((*curspecies)->organisms).begin()))->orig_fitness)>
		highest_fitness1) {
			highest_fitness1=((*(((*curspecies)->organisms).begin()))->orig_fitness);
			highest_last_changed1=0;
			std::cout<<"NEW POPULATION RECORD FITNESS: "<<highest_fitness1<<std::endl;
		}
	else {
		++highest_last_changed1;
		std::cout<<highest_last_changed1<<" generations since last population fitness record: "<<highest_fitness1<<std::endl;
	}


	//Check for stagnation- if there is stagnation, perform delta-coding
	if (highest_last_changed1>=NEAT::dropoff_age+5) {

		//    cout<<"PERFORMING DELTA CODING"<<endl;

                std::cout<<"Aditya: Delta coding in population.cpp->epoch"<<std::endl; //Aditya:for debug
		highest_last_changed1=0;

		half_pop=NEAT::pop_size/2;

		//    cout<<"half_pop"<<half_pop<<" pop_size-halfpop: "<<pop_size-half_pop<<endl;

		curspecies=sorted_species.begin();

		(*(((*curspecies)->organisms).begin()))->super_champ_offspring=half_pop;
		(*curspecies)->expected_offspring=half_pop;
		(*curspecies)->age_of_last_improvement=(*curspecies)->age;

		++curspecies;

		if (curspecies!=sorted_species.end()) {

			(*(((*curspecies)->organisms).begin()))->super_champ_offspring=NEAT::pop_size-half_pop;
			(*curspecies)->expected_offspring=NEAT::pop_size-half_pop;
			(*curspecies)->age_of_last_improvement=(*curspecies)->age;

			++curspecies;

			//Get rid of all species under the first 2
			while(curspecies!=sorted_species.end()) {
				(*curspecies)->expected_offspring=0;
				++curspecies;
			}
		}
		else {
			curspecies=sorted_species.begin();
			(*(((*curspecies)->organisms).begin()))->super_champ_offspring+=NEAT::pop_size-half_pop;
			(*curspecies)->expected_offspring=NEAT::pop_size-half_pop;
		}

	}
	//STOLEN BABIES:  The system can take expected offspring away from
	//  worse species and give them to superior species depending on
	//  the system parameter babies_stolen (when babies_stolen > 0)
	else if (NEAT::babies_stolen>0) {
		//Take away a constant number of expected offspring from the worst few species

		stolen_babies=0;
		curspecies=sorted_species.end();
		curspecies--;
		while ((stolen_babies<NUM_STOLEN)&&
			(curspecies!=sorted_species.begin())) {

				//cout<<"Considering Species "<<(*curspecies)->id<<": age "<<(((*curspecies)->age))<<" expected offspring "<<(((*curspecies)->expected_offspring))<<endl;

				if ((((*curspecies)->age)>5)&&
					(((*curspecies)->expected_offspring)>2)) {
						//cout<<"STEALING!"<<endl;

						//This species has enough to finish off the stolen pool
						if (((*curspecies)->expected_offspring-1)>=(NUM_STOLEN-stolen_babies)) {
							(*curspecies)->expected_offspring-=(NUM_STOLEN-stolen_babies);
							stolen_babies=NUM_STOLEN;
						}
						//Not enough here to complete the pool of stolen
						else {
							stolen_babies+=(*curspecies)->expected_offspring-1;
							(*curspecies)->expected_offspring=1;

						}
					}

					curspecies--;

					//if (stolen_babies>0)
					//cout<<"stolen babies so far: "<<stolen_babies<<endl;
			}

			//cout<<"STOLEN BABIES: "<<stolen_babies<<endl;

			//Mark the best champions of the top species to be the super champs
			//who will take on the extra offspring for cloning or mutant cloning
			curspecies=sorted_species.begin();

			//Determine the exact number that will be given to the top three
			//They get , in order, 1/5 1/5 and 1/10 of the stolen babies
			one_fifth_stolen=NEAT::babies_stolen/5;
			one_tenth_stolen=NEAT::babies_stolen/10;

			//Don't give to dying species even if they are champs
			while((curspecies!=sorted_species.end())&&((*curspecies)->last_improved()>NEAT::dropoff_age))
				++curspecies;

			//Concentrate A LOT on the number one species
			if ((stolen_babies>=one_fifth_stolen)&&(curspecies!=sorted_species.end())) {
				(*(((*curspecies)->organisms).begin()))->super_champ_offspring=one_fifth_stolen;
				(*curspecies)->expected_offspring+=one_fifth_stolen;
				stolen_babies-=one_fifth_stolen;
				//cout<<"Gave "<<one_fifth_stolen<<" babies to Species "<<(*curspecies)->id<<endl;
				//      cout<<"The best superchamp is "<<(*(((*curspecies)->organisms).begin()))->gnome->genome_id<<endl;

				//Print this champ to file "champ" for observation if desired
				//IMPORTANT:  This causes generational file output 
				//print_Genome_tofile((*(((*curspecies)->organisms).begin()))->gnome,"champ");

				curspecies++;

			}

			//Don't give to dying species even if they are champs
			while((curspecies!=sorted_species.end())&&((*curspecies)->last_improved()>NEAT::dropoff_age))
				++curspecies;

			if ((curspecies!=sorted_species.end())) {
				if (stolen_babies>=one_fifth_stolen) {
					(*(((*curspecies)->organisms).begin()))->super_champ_offspring=one_fifth_stolen;
					(*curspecies)->expected_offspring+=one_fifth_stolen;
					stolen_babies-=one_fifth_stolen;
					//cout<<"Gave "<<one_fifth_stolen<<" babies to Species "<<(*curspecies)->id<<endl;
					curspecies++;

				}
			}

			//Don't give to dying species even if they are champs
			while((curspecies!=sorted_species.end())&&((*curspecies)->last_improved()>NEAT::dropoff_age))
				++curspecies;

			if (curspecies!=sorted_species.end())
				if (stolen_babies>=one_tenth_stolen) {
					(*(((*curspecies)->organisms).begin()))->super_champ_offspring=one_tenth_stolen;
					(*curspecies)->expected_offspring+=one_tenth_stolen;
					stolen_babies-=one_tenth_stolen;

					//cout<<"Gave "<<one_tenth_stolen<<" babies to Species "<<(*curspecies)->id<<endl;
					curspecies++;

				}

				//Don't give to dying species even if they are champs
				while((curspecies!=sorted_species.end())&&((*curspecies)->last_improved()>NEAT::dropoff_age))
					++curspecies;

				while((stolen_babies>0)&&
					(curspecies!=sorted_species.end())) {
						//Randomize a little which species get boosted by a super champ

						if (randfloat()>0.1)
							if (stolen_babies>3) {
								(*(((*curspecies)->organisms).begin()))->super_champ_offspring=3;
								(*curspecies)->expected_offspring+=3;
								stolen_babies-=3;
								//cout<<"Gave 3 babies to Species "<<(*curspecies)->id<<endl;
							}
							else {
								//cout<<"3 or less babies available"<<endl;
								(*(((*curspecies)->organisms).begin()))->super_champ_offspring=stolen_babies;
								(*curspecies)->expected_offspring+=stolen_babies;
								//cout<<"Gave "<<stolen_babies<<" babies to Species "<<(*curspecies)->id<<endl;
								stolen_babies=0;

							}

							curspecies++;

							//Don't give to dying species even if they are champs
							while((curspecies!=sorted_species.end())&&((*curspecies)->last_improved()>NEAT::dropoff_age))
								++curspecies;

					}

					//cout<<"Done giving back babies"<<endl;

					//If any stolen babies aren't taken, give them to species #1's champ
					if (stolen_babies>0) {

						//cout<<"Not all given back, giving to best Species"<<endl;

						curspecies=sorted_species.begin();
						(*(((*curspecies)->organisms).begin()))->super_champ_offspring+=stolen_babies;
						(*curspecies)->expected_offspring+=stolen_babies;
						stolen_babies=0;
					}

	}


	//Kill off all Organisms marked for death.  The remainder
	//will be allowed to reproduce.
        int count_org = 0;
	curorg=organisms.begin();
	while(curorg!=organisms.end()) {
                count_org = count_org + 1;
		if (((*curorg)->eliminate)) {
                        std::cout << "Aditya: Error cannot eliminate organisms"<<std::endl;


			//Remove the organism from its Species
			((*curorg)->species)->remove_org(*curorg);

			//Delete the organism from memory
			delete (*curorg);

			//Remember where we are
			deadorg=curorg;
			++curorg;

			//iter2 =  v.erase(iter); 

			//Remove the organism from the master list
			curorg=organisms.erase(deadorg);

		}
		else {
			++curorg;
		}

	}
        std::cout << "Aditya: Number of organisms"<<count_org<<std::endl;

	//cout<<"Reproducing"<<endl;

	//Perform reproduction.  Reproduction is done on a per-Species
	//basis.  (So this could be paralellized potentially.)
	//	for(curspecies=species.begin();curspecies!=species.end();++curspecies) {

	//KENHACK                                                                      
	//		for(std::vector<Species*>::iterator curspecies2=species.begin();curspecies2!=species.end();++curspecies2) {
	//		  std::cout<<"PRE in repro specloop SPEC EXISTING number "<<(*curspecies2)->id<<std::endl;
	//	}

	//	(*curspecies)->reproduce(generation,this,sorted_species);


	//}    


	curspecies=species.begin();
	int last_id=(*curspecies)->id;
	while(curspecies!=species.end()) {
	  (*curspecies)->reproduce(generation,this,sorted_species);

	  //Set the current species to the id of the last species checked
	  //(the iterator must be reset because there were possibly vector insertions during reproduce)
	  std::vector<Species*>::iterator curspecies2=species.begin();
	  while(curspecies2!=species.end()) {
	    if (((*curspecies2)->id)==last_id)
	      curspecies=curspecies2;
	    curspecies2++;
	  }

	  //Move to the next on the list
	  curspecies++;
	  
	  //Record where we are
	  if (curspecies!=species.end())
	    last_id=(*curspecies)->id;

	}

	//cout<<"Reproduction Complete"<<endl;


	//Destroy and remove the old generation from the organisms and species  
	curorg=organisms.begin();
	while(curorg!=organisms.end()) {

	  //Remove the organism from its Species
	  ((*curorg)->species)->remove_org(*curorg);

	  //std::cout<<"deleting org # "<<(*curorg)->gnome->genome_id<<std::endl;

	  //Delete the organism from memory
	  delete (*curorg);
	  
	  //Remember where we are
	  deadorg=curorg;
	  ++curorg;
	  
	  //std::cout<<"next org #  "<<(*curorg)->gnome->genome_id<<std::endl;

	  //Remove the organism from the master list
	  curorg=organisms.erase(deadorg);

	  //std::cout<<"nnext org # "<<(*curorg)->gnome->genome_id<<std::endl;

	}

	//Remove all empty Species and age ones that survive
	//As this happens, create master organism list for the new generation
	curspecies=species.begin();
	orgcount=0;
	while(curspecies!=species.end()) {
		if (((*curspecies)->organisms.size())==0) {
			delete (*curspecies);

			deadspecies=curspecies;
			++curspecies;

			curspecies=species.erase(deadspecies);
		}
		//Age surviving Species and 
		//Rebuild master Organism list: NUMBER THEM as they are added to the list
		else {
			//Age any Species that is not newly created in this generation
			if ((*curspecies)->novel) {
				(*curspecies)->novel=false;
			}
			else ++((*curspecies)->age);

			//Go through the organisms of the curspecies and add them to 
			//the master list
			for(curorg=((*curspecies)->organisms).begin();curorg!=((*curspecies)->organisms).end();++curorg) {
				((*curorg)->gnome)->genome_id=orgcount++;
				organisms.push_back(*curorg);
			}
			++curspecies;

		}
	}      

	//Remove the innovations of the current generation
	curinnov=innovations.begin();
	while(curinnov!=innovations.end()) {
		delete (*curinnov);

		deadinnov=curinnov;
		++curinnov;

		curinnov=innovations.erase(deadinnov);
	}

	//DEBUG: Check to see if the best species died somehow
	// We don't want this to happen
	curspecies=species.begin();
	best_ok=false;
	while(curspecies!=species.end()) {
		if (((*curspecies)->id)==best_species_num) best_ok=true;
		++curspecies;
	}
	if (!best_ok) {
		//cout<<"ERROR: THE BEST SPECIES DIED!"<<endl;
	}
	else {
		//cout<<"The best survived: "<<best_species_num<<endl;
	}

	//DEBUG: Checking the top organism's duplicate in the next gen
	//This prints the champ's child to the screen
	for(curorg=organisms.begin();curorg!=organisms.end();++curorg) {
		if ((*curorg)->pop_champ_child) {
			//cout<<"At end of reproduction cycle, the child of the pop champ is: "<<(*curorg)->gnome<<endl;
		}
	}

	//cout<<"babies_stolen at end: "<<babies_stolen<<endl;

	//cout<<"Epoch complete"<<endl; 

	return true;

}

bool Population::rank_within_species() {
	std::vector<Species*>::iterator curspecies;

	//Add each Species in this generation to the snapshot
	for(curspecies=species.begin();curspecies!=species.end();++curspecies) {
		(*curspecies)->rank();
	}

	return true;
}
