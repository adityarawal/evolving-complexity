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
#include "organism.h"
#include <iostream>

using namespace NEAT;

Organism::Organism(double fit1, double fit2, Genome *g,int gen, const char* md) {
	fitness1=fit1;//(Aditya: for NSGA-2)
	fitness2=fit2;//(Aditya: for NSGA-2)
	orig_fitness=fitness1;//(Aditya: for NSGA-2)
        front_num = -1;//(Aditya: for NSGA-2)
        crowd_dist = 0.0;//(Aditya: for NSGA-2)
	evaluated = false; //Aditya: for speed-up by preventing re-evaluation of the elites
        gnome=g;
	net=gnome->genesis(gnome->genome_id);
	species=0;  //Start it in no Species
	expected_offspring=0;
	generation=gen;
	eliminate=false;
	error=0;
	winner=false;
	champion=false;
	super_champ_offspring=0;

	// If md is null, then we don't have metadata, otherwise we do have metadata so copy it over
	if(md == 0) {
		strcpy(metadata, "");
	} else {
		strncpy(metadata, md, 128);
	}

	time_alive=0;

	//DEBUG vars
	pop_champ=false;
	pop_champ_child=false;
	high_fit=0;
	mut_struct_baby=0;
	mate_baby=0;

	modified = true;
}

Organism::Organism(const Organism& org)
{
	fitness1 = org.fitness1; //(Aditya: for NSGA-2)
	fitness2 = org.fitness2;//(Aditya: for NSGA-2)
	orig_fitness = org.orig_fitness;//(Aditya: for NSGA-2)
        front_num = org.front_num; //(Aditya: for NSGA-2)
        crowd_dist = org.crowd_dist; //(Aditya: for NSGA-2)
	evaluated = false; //Aditya: for speed-up by preventing re-evaluation of the elites
	gnome = new Genome(*(org.gnome));	// Associative relationship
	//gnome = org.gnome->duplicate(org.gnome->genome_id);
	net = new Network(*(org.net)); // Associative relationship
	species = org.species;	// Delegation relationship
	expected_offspring = org.expected_offspring;
	generation = org.generation;
	eliminate = org.eliminate;
	error = org.error;
	winner = org.winner;
	champion = org.champion;
	super_champ_offspring = org.super_champ_offspring;

	strcpy(metadata, org.metadata);
	//printf("copying %s did it work? %s", org.metadata, metadata);

	time_alive = org.time_alive;
	pop_champ = org.pop_champ;
	pop_champ_child = org.pop_champ_child;
	high_fit = org.high_fit;
	mut_struct_baby = org.mut_struct_baby;
	mate_baby = org.mate_baby;

	modified = false;
}

Organism::~Organism() {
	delete net;
	delete gnome;
}

void Organism::remove_inactive_genome() {
        std::vector<Gene*>::iterator curgene;
        std::vector<NNode*>::iterator curnode;
        std::vector<Link*>::iterator curlink;
        std::vector<int> active_node_ids;
	std::vector<int>::iterator location;
	std::vector<NNode*>::iterator deadnode;
	std::vector<Gene*>::iterator deadgene;

        //Recursively find all the valid paths between input and output
        net->find_active_paths();

        //Create a list of all active node_ids
        //This is required because genome->nodes are 
        //not the same as net->all_nodes
        //Also, find and mark the active links in the genome
	for(curnode=net->all_nodes.begin();curnode!=net->all_nodes.end();++curnode) {
                if ((*curnode)->on_active_path) {
                        active_node_ids.push_back((*curnode)->node_id);
                }
                for(curlink=((*curnode)->incoming).begin();curlink!=((*curnode)->incoming).end();++curlink) {
                        if ((*curlink)->on_active_path) {
                                //Find the corresponding link in the genome
                                //link and mark its on_active_path
                                for(curgene=gnome->genes.begin();curgene!=gnome->genes.end();++curgene) {
                                        if (((*curgene)->lnk->in_node->node_id == (*curlink)->in_node->node_id) &&
                                            ((*curgene)->lnk->out_node->node_id == (*curlink)->out_node->node_id) &&
                                            ((*curgene)->lnk->is_recurrent == (*curlink)->is_recurrent) &&
                                            ((*curgene)->lnk->link_gtype == (*curlink)->link_gtype)) {
                                                  (*curgene)->lnk->on_active_path = true;
                                                  break;
                                        }
                                }

                        }
                }
        }
        //Delete inactive nodes
     	curnode=gnome->nodes.begin();
	while(curnode!=gnome->nodes.end()) {
                location = std::find(active_node_ids.begin(),active_node_ids.end(),(*curnode)->node_id);
                if (location==active_node_ids.end()){
                        delete (*curnode);
	                //Remember where we are
	                deadnode=curnode;
	                ++curnode;
			curnode=gnome->nodes.erase(deadnode);
                }
		else {
			++curnode;
		}

	}    

        //Delete inactive genes
     	curgene=gnome->genes.begin();
	while(curgene!=gnome->genes.end()) {
                if (!((*curgene)->lnk->on_active_path)){
                        delete (*curgene);
	                //Remember where we are
	                deadgene=curgene;
	                ++curgene;
			curgene=gnome->genes.erase(deadgene);
                }
		else {
			++curgene;
		}

	}    
	
}

void Organism::update_phenotype() {

	//First, delete the old phenotype (net)
	delete net;

	//Now, recreate the phenotype off the new genotype
	net=gnome->genesis(gnome->genome_id);

	modified = true;
}

bool Organism::print_to_file(char *filename) {
	
	std::ofstream oFile(filename);

	return write_to_file(oFile);
}

bool Organism::write_to_file(std::ostream &outFile) {
	
	char tempbuf2[1024];
	if(modified == true) {
		sprintf(tempbuf2, "/* Organism #%d Fitness: %f Time: %d */\n", (gnome)->genome_id, fitness1, fitness2, time_alive);
	} else {
		sprintf(tempbuf2, "/* %s */\n", metadata);
	}
	outFile << tempbuf2;
	gnome->print_to_file(outFile);
	return 1;
}

//// Print the Organism's genome to a file preceded by a comment             
////   detailing the organism's species, number, and fitness
//bool Organism::print_to_file(char *filename) {
//
//ofstream oFile(filename,ios::out);
//
//cout<<"FILENAME: "<<filename<<endl;
//
////Make sure it worked
//if (!oFile) {
//cerr<<"Can't open "<<filename<<" for output"<<endl;
//return 0;
//}
//
////Put the fitness and other information for each organism in a comment
//oFile<<endl<<"/* Organism #"<<gnome->genome_id<<" Fitness: "<<fitness<<" *///"<<endl;
//
//gnome->print_to_file(oFile);
//
//return 1;
//}

bool NEAT::order_orgs(Organism *x, Organism *y) {//(Aditya: for NSGA-2)
	return (x)->fitness1 > (y)->fitness1;
}

bool NEAT::order_orgs_by_adjusted_fit(Organism *x, Organism *y) {//(Aditya: for NSGA-2)
	return (x)->fitness1 / (x->species)->organisms.size()  > (y)->fitness1 / (y->species)->organisms.size();
}
