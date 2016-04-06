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
#include "experiments.h"
#include <cstring>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <iterator>     // std::ostream_iterator
#include <algorithm>    // std::copy

int toBin(double v, int num_bin) {
	if(v>0.999) {
                v=0.999;
        }
        else {
                v =v+0.0001;//To avoid spurious error due to inaccuracy of double
        }
        return (int)((v)*(double)num_bin);//type casting to integer floors the number
}

void ConvertToBinary(int n, std::vector <int> &sequence)
{
    if (n / 2 != 0) {
        ConvertToBinary(n / 2, sequence);
    }
    sequence.push_back(n%2);
}

void read_input_file_data(std::vector < vector < double > > &input_data, int num_input_nodes){

    cout<<"Reading Input Sequence Data from file "<<endl;
    ifstream XtrainFile("input_data.txt",ios::in);//File name("/scratch/cluster/aditya/DBN_research/rp_deep/original_code/F.txt",ios::in)
    double d;
    std::string lineData;

    while (getline(XtrainFile, lineData)) {
            std::vector < double > input_row;
            std::stringstream lineStream(lineData);
            int feature_count = 1;

            input_row.push_back(1.0); //Inserting 1.0 for Bias node
            while ((lineStream >> d) && (feature_count <= num_input_nodes)) {//Read only as many features as there are input nodes
                    feature_count += 1;
                    input_row.push_back(d);
            }
            input_data.push_back(input_row);
    }
    XtrainFile.close();
    cout<<"Completed Reading Sequence Data, Number of sequences: "<<input_data.size()<<", Number of Features in each Image (Including Bias input): "<<input_data[0].size()<<endl;
    cout<<"Also Inserted Ones at the start of each input image feature for bias node "<<endl;
}

void generate_real_input_data(std::vector < vector < double > > &input_data, int num_input_nodes, int num_sequences){

    cout<<"Generating Real Input Sequence Data "<<endl;
    input_data.resize(num_sequences, std::vector<double>());
    for (int i=0; i <num_sequences; i++) {
            for (int j=0; j<NEAT::input_sequence_len;j++) {
                    double rand_num = randfloat();
                    input_data[i].push_back(rand_num);
            }
    }
}

void generate_bin_input_data(std::vector < vector < double > > &input_data, int num_input_nodes, int num_sequences, int num_bin){

    cout<<"Generating Binned Input Sequence Data (range (1.0/num_bin) - 1.0) with "<< num_bin <<" bins"<<endl;
    input_data.resize(num_sequences, std::vector<double>());
    for (int i=0; i <num_sequences; i++) {
            for (int j=0; j<NEAT::input_sequence_len;j++) {
                    double rand_num = randfloat();
                    int bin_num = toBin(rand_num, num_bin)+1; //+1 makes sure that the value ranges from 1-num+bin
                    double scaled_bin_num = bin_num/num_bin; //making range (1.0/num_bin) - 1.0
                    input_data[i].push_back(scaled_bin_num);
            }
    }
}

void generate_discrete_input_data(std::vector < vector < double > > &input_data){
    
    cout<<"Generating Discrete Input Sequence Data "<<endl;
    int num_sequences = pow(2, NEAT::input_sequence_len);
    input_data.resize(num_sequences, std::vector<double>());

    std::vector <int> sequence;
    for (int i=0; i <num_sequences; i++) {
            sequence.clear();
            ConvertToBinary(i, sequence);
            if (sequence.size()<NEAT::input_sequence_len){//If sequence length is less than the required
                    for (int j=0; j<NEAT::input_sequence_len-sequence.size();j++) {
                            input_data[i].push_back(0);
                    }
            }
            for (int j=0; j < sequence.size(); j++) {
                    input_data[i].push_back(sequence[j]);
            }
    }
    for (int i=0; i<input_data.size(); i++) {
            for (int j=0; j<input_data[i].size(); j++) {
                    if(input_data[i][j] == 0) {//Swap 0s with -1s
                            input_data[i][j] = -1;
                    }
            }
    }
} 

void generate_Tmaze_input_data(int num_input_nodes, std::vector < vector < double > > &input_data, int num_trials, int aisle_length){
    cout<<"Generating Tmaze Input Data "<<endl;
    int trial_length = aisle_length*NEAT::input_sequence_len;
    input_data.resize(num_trials*trial_length, std::vector<double>());
    std::vector < vector < double > > input_sequences; //All the input sequences (does not include other inputs like wall distance, valid etc.

    //Generate all possible binary (1/-1) sequences of a given length 
    generate_discrete_input_data(input_sequences);
    
    for (int i = 0; i<num_trials; i++) {
            
            //At the start of the trial, input sequence is shown to the agent  
            bool onetime = false;//Used to make sure that light sequence input is given only once at the start of the trial
            for (int j=0; j<NEAT::input_sequence_len; j++) {
                    int Tjunction_dist = aisle_length;
                    for (int k=0; k<aisle_length; k++ ) {
                            Tjunction_dist = Tjunction_dist - 1;
                            input_data[i*trial_length+j*aisle_length+k].push_back(Tjunction_dist);
                            if (!onetime) {//If the light have not been shown to the agent (Basically, start of the trial)
                                    if (k<input_sequence_len) {
                                            double rand_num = randfloat();
                                            //if (i%2 == 0.0) {
                                            input_data[i*trial_length+j*aisle_length+k].push_back(input_sequences[i][k]);
                                            //}
                                            //else {
                                            //}
                                    }
                                    else {
                                            input_data[i*trial_length+j*aisle_length+k].push_back(0);

                                    }
                            }
                            else { //No more light sequence input on subsequent turns
                                    input_data[i*trial_length+j*aisle_length+k].push_back(0);
                                    
                            }
                    }
                    onetime = true; //Set it so that the lights are not shown again (in subsequent T-junctions)
            }
    }
}

void freeze_update_genome(Genome *start_genome) {
         std::cout<<"Freezing start genome"<<std::endl;

         //Freeze the winning genome so that new stuff can be added to it
         start_genome->freeze_genome(); //Freeze current genome to prevent any new incoming connections to the existing nodes and any weight changes on this part of the network
         
         //Find the last innovation number
         double max_innov_num =0.0;
         std::vector<Gene*>::iterator curgene;
         for(curgene=start_genome->genes.begin();curgene!=start_genome->genes.end();++curgene) {
         	if ((*curgene)->innovation_num > max_innov_num) {
                     max_innov_num = (*curgene)->innovation_num;
             }
         }
         max_innov_num = max_innov_num + 1.0;//points to the newly added node

         //Add #NEAT::batch_size new output nodes
         std::cout<<"Adding "<<batch_size<<" new output nodes to the frozen genome"<<std::endl;
         start_genome->add_output_nodes(NEAT::batch_size, max_innov_num);
}

//#define NO_SCREEN_OUT 
//Perform evolution on digit recognition task, for gens generations
Population *memory_test(int gens) {
    Population *pop=0;
    Genome *start_genome;
    Genome *new_winner_genome = NULL; //Current winner 
    Genome *last_winner_genome = NULL; //Last winner with lesser output nodes
    char curword[20];
    int id;

    ostringstream *fnamebuf;
    int gen;
 
    int evals[NEAT::num_runs];  //Hold records for each run
    int genes[NEAT::num_runs];
    int nodes[NEAT::num_runs];
    int winnernum;
    int winnergenes;
    int winnernodes;
    //For averaging
    int totalevals=0;
    int totalgenes=0;
    int totalnodes=0;
    int expcount;
    int samples;  //For averaging
    int before_nw_size, after_nw_size;

    int current_output_nodes = 0;
    vector<Organism*>::iterator curorg;
    
    
    memset (evals, 0, NEAT::num_runs * sizeof(int));
    memset (genes, 0, NEAT::num_runs * sizeof(int));
    memset (nodes, 0, NEAT::num_runs * sizeof(int));

      
    //Run experiments
    for(expcount=0;expcount<NEAT::num_runs;expcount++) {

          std::vector < vector < double > > input_data; //Image data
          vector< vector <double> > independent_archive;//Archive of the independent output features (num_frozen_output_nodes x input_data.size())

          ifstream iFile("memory_startgenes",ios::in);
          cout<<"START Memory Evolution TEST"<<endl;
          cout<<"Reading in the start genome"<<endl;
          //Read in the start Genome
          iFile>>curword;
          iFile>>id;
          cout<<"Reading in Genome id "<<id<<endl;
          start_genome=new Genome(id,iFile);
          iFile.close();

          //Count number of input nodes, so that only so many features are read from the file
          int num_input_nodes = 0;
          for (int i=0; i<start_genome->nodes.size(); i++) {
                  if (start_genome->nodes[i]->gen_node_label==1) {//Input node
                          num_input_nodes += 1;
                  }
                  if (start_genome->nodes[i]->gen_node_label==2) {//Output node
                          current_output_nodes += 1;
                  }
          }
          std::cout << "Number of Start Genome Input nodes: "<< num_input_nodes<<std::endl;    
          std::cout << "Number of Start Genome Output nodes: "<< current_output_nodes<<std::endl;    

          //Freeze the startgenome if required
          if (NEAT::frozen_startgenome==1) {
               freeze_update_genome(start_genome);
               //start_genome->freeze_genome(); //Freeze current genome to prevent any new incoming connections to the existing nodes and any weight changes on this part of the network
          }
          //Resetting the pointers so that each run/experiment has new winners
          new_winner_genome = NULL;
          last_winner_genome = NULL;
          
          //Clear the archive for a new run
          independent_archive.clear();

          //Read the Image data and output labels (digits) from text file
          int num_trials = pow(2, NEAT::input_sequence_len); //Number of Trials
          int aisle_length = 10; //Aisle length (Fixed for now)
          int max_rand_activate = 20; //Not used currently
          generate_Tmaze_input_data(num_input_nodes, input_data, num_trials, aisle_length);
          int max_nw_size = 7;
          //Spawn the Population
          cout<<"Spawning Population off Genome "<<start_genome->genome_id<<endl;;

          pop=new Population(start_genome,NEAT::pop_size);
          
          cout<<"Verifying Spawned Pop"<<endl;
          pop->verify();

          for (gen=1;gen<=gens;gen++) {

                  cout<<"Epoch "<<gen<<endl;	

                  current_output_nodes = pop->organisms[0]->net->outputs.size();//Current output nodes in each member of the population

                  char temp[50];
                  sprintf (temp, "gen_%d_%d", gen,current_output_nodes);

                  bool success = false;
                  success = memory_epoch(pop,gen,temp,winnernum,winnergenes,winnernodes, input_data, independent_archive, num_trials, aisle_length, max_nw_size);
                  //Check for success
                  if  (success && (current_output_nodes == NEAT::max_output_nodes)) {//If all independent output features have been discovered
                          //Collect Stats on end of experiment
                          evals[expcount]=NEAT::pop_size*(gen-1)+winnernum;
                          genes[expcount]=winnergenes;
                          nodes[expcount]=winnernodes;
                          gen=gens;
                          delete last_winner_genome;
                          delete new_winner_genome;
                  }
                  else if (success && (current_output_nodes < NEAT::max_output_nodes)) { //If current outputs (<total outputs) are independent)

                          //Grab the last winning genome from the population and make sure its frozen part is unchanged
                          int winner_genomeid;
                          bool check_genome;
                          for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
                            if ((*curorg)->winner) {
                              before_nw_size = (((*curorg)->gnome)->compute_genome_size());
                              (*curorg)->remove_inactive_genome(); //Delete inactive parts (nodes/genes) of the genome
                              after_nw_size = (((*curorg)->gnome)->compute_genome_size());
                              if (before_nw_size!=after_nw_size) {
                                      std::cout<<"ERROR: in remove_inactive_genome()"<<std::endl;
                                      exit(0);
                              }
                              winner_genomeid = ((*curorg)->gnome)->genome_id;
                              new_winner_genome=((*curorg)->gnome)->duplicate(1);
                              //break;
                            }
                          }

                          //With each newly added LSTM output,
                          //the maximum network size needs to be updated 
                          max_nw_size += 4; 


                          //Verify that the frozen part of the network is unchanged
                          if (last_winner_genome != NULL) {//First time winning genome is created and frozen
                                  check_genome = last_winner_genome->compare_frozen_genome(new_winner_genome);
                                  if (!check_genome) {
                                          std::cout<<"Frozen genome has been modified, exiting.... "<<std::endl;
                                          exit(0);
                                  }
                          }

                          //Save the genome for comparison with the next winner
                          last_winner_genome=new_winner_genome->duplicate(1);

                          //Freeze the winning genome so that new stuff can be added to it
                          new_winner_genome->freeze_genome(); //Freeze current genome to prevent any new incoming connections to the existing nodes and any weight changes on this part of the network

                          //Add #NEAT::batch_size new output nodes
                          new_winner_genome->add_output_nodes(NEAT::batch_size, pop->cur_innov_num);


                          delete pop;//Is this required??
                          
                          //RE-SPAWN THE POPULATION
                          cout<<"Spawning Population off Genome "<<winner_genomeid<<" with "<<current_output_nodes <<" frozen + "<<NEAT::batch_size<<" new output nodes"<<endl;
                          pop=new Population(new_winner_genome,NEAT::pop_size);
                          cout<<"Verifying Spawned Pop"<<endl;
                          pop->verify();
                          //pop->print_to_file_by_species(strcat(temp, "_del"));
                  }
                  
          }

          if (expcount<NEAT::num_runs-1) delete pop;
      
    }

    //Average and print stats
    cout<<"Nodes: "<<endl;
    for(expcount=0;expcount<NEAT::num_runs;expcount++) {
      cout<<nodes[expcount]<<endl;
      totalnodes+=nodes[expcount];
    }
    
    cout<<"Genes: "<<endl;
    for(expcount=0;expcount<NEAT::num_runs;expcount++) {
      cout<<genes[expcount]<<endl;
      totalgenes+=genes[expcount];
    }
    
    cout<<"Evals "<<endl;
    samples=0;
    for(expcount=0;expcount<NEAT::num_runs;expcount++) {
      cout<<evals[expcount]<<endl;
      if (evals[expcount]>0)
      {
	totalevals+=evals[expcount];
	samples++;
      }
    }

    cout<<"Failures: "<<(NEAT::num_runs-samples)<<" out of "<<NEAT::num_runs<<" runs"<<endl;
    cout<<"Average Nodes: "<<(samples>0 ? (double)totalnodes/samples : 0)<<endl;
    cout<<"Average Genes: "<<(samples>0 ? (double)totalgenes/samples : 0)<<endl;
    cout<<"Average Evals: "<<(samples>0 ? (double)totalevals/samples : 0)<<endl;

    return pop;

}
double BintoValue(int v, int num_bin) {
        return (double)(((double)(v)/(double)num_bin)+0.0001);//type casting to integer floors the number
}

double log_anybase(double value, int base) {
        return (double) log(value)/log((double) base);
}

class histogram
{
        std::vector<int> bin_count;
        std::vector<int> num_bins;
        std::vector<double> granularity;
        int total;
  
        public:
        void initialize1d(std::vector<double> x, int x_num_bins) {
                num_bins.push_back(x_num_bins);
        	granularity.push_back((double)(1/num_bins[0]));
        	for (int i=0; i<num_bins[0]; i++) { 
                        bin_count.push_back(0);  //Initializing
                }
                for (int i=0; i<x.size(); i++) {
                        bin_count[toBin(x[i], num_bins[0])]++;//Incrementing bin count
                }
        	total=x.size();
        }
        void initialize2d(std::vector<double> x, std::vector<double> y,int x_num_bins, int y_num_bins) {
        	num_bins.push_back(x_num_bins);
        	num_bins.push_back(y_num_bins);
                granularity.push_back(1.0/(double)x_num_bins);
                granularity.push_back(1.0/(double)y_num_bins);
        	for (int i=0; i<num_bins[0]*num_bins[1]; i++) { 
                        bin_count.push_back(0);  //Initializing
                }
        	for(int i=0;i<x.size();i++) {
        		bin_count[toBin2d(x[i],y[i])]++;
        	}
        	total=x.size();
        }
        int toBin2d(double x,double y) {
        	if(x>0.999) x=0.999;
        	if(y>0.999) y=0.999;
        	return toBin(x, num_bins[0])*num_bins[1]+toBin(y, num_bins[1]);
        }
        
        double probability(double x,double y) {
        	return ((double)bin_count[toBin2d(x,y)])/((double)total);///(granularity[0]*granularity[1]);
        }
        
        //int toBin(double v, int num_bin) {
        //	if(v==1.0) {
        //                v=0.999;
        //        }
        //        else {
        //                v =v+0.0001;//To avoid spurious error due to inaccuracy of double
        //        }
        //        return (int)((v)*(double)num_bin);//type casting to integer floors the number
        //}
        double probability(double v) {
        	return ((double)bin_count[toBin(v, num_bins[0])])/((double)total);///(granularity[0]);
        }
        void print_bin_counts() {
                if (num_bins.size() > 1) {
                        for (int i=0; i<num_bins[0]*num_bins[1]; i++) {
                                std::cout<<i<<" "<<bin_count[i]<<std::endl;
                        }
                        std::cout<<endl;
                }
                else {
                        for (int i=0; i<num_bins[0]; i++) {
                                std::cout<<i<<" "<<bin_count[i]<<std::endl;
                        }
                        std::cout<<endl;
                }
        }
        
        static double standard_deviation(std::vector<double> xvals, int x_num_bins) {

        	histogram hx;// = new histogram();
        	hx.initialize1d(xvals,x_num_bins);
                //hx.print_bin_counts();
                double px; //Probability
                double mean = 0.0;
                double variance = 0.0;
                double std_dev = 0.0;
                double granularity_x = (double) (1.0/(double) x_num_bins); 

        	for(double x=0.0000001;x<1.0;x=x+granularity_x) {//Starting with x=0.0000001 to avoid double inaccuracy
                        px=hx.probability(x);
                        mean += px*x;
                }
        	for(double x=0.0000001;x<1.0;x=x+granularity_x) {//Starting with x=0.0000001 to avoid double inaccuracy
                        px=hx.probability(x);
                        variance += px*(pow((x-mean),2));
                }
                if (variance > 0.0) {
                        std_dev = sqrt(variance);
                        //std::cout<<" Standard Deviation: "<<std_dev<<" "<<sqrt(0.021)<<endl;
                }
                else {
                        std_dev = 0.0;
                }
                return std_dev;
        }
        
        static double entropy(std::vector<double> xvals, int x_num_bins) {
        	histogram hx;// = new histogram();
        	hx.initialize1d(xvals,x_num_bins);
                //hx.print_bin_counts();
                double granularity_x = (double) (1.0/(double) x_num_bins); 
                double entropy=0.0;
        	for(double x=0.0000001;x<1.0;x=x+granularity_x) {//Starting with x=0.0000001 to avoid double inaccuracy
                        double px=hx.probability(x);
                        if (px == 1) {//Diversity experiment (to ensure neuron value is not fixed)
                                entropy = 0.0;
                        }
                        else if(px!=0.0) {
        			entropy+=(-1.0)*px*(log_anybase(px, x_num_bins));
                        }
                }
                //std::cout<<"Entropy:"<<entropy<<std::endl;
                if (entropy > 1.0) {
                        std::cout<<" Entropy IS MORE THAN : 1.0 "<<entropy<<std::endl;//Uniform distb shd have max = 1.0
                        exit(0);
                }
        	return entropy;		

        }
        static double mutual_inf(std::vector<double> xvals, std::vector<double> yvals, int x_num_bins, int y_num_bins) {
        	histogram hxy;// = new histogram();
        	histogram hx;// = new histogram();
        	histogram hy;// = new histogram();
        	hx.initialize1d(xvals,x_num_bins);
        	hy.initialize1d(yvals,y_num_bins);
        	hxy.initialize2d(xvals,yvals,x_num_bins, y_num_bins);
                
                //hx.print_bin_counts();
                //hy.print_bin_counts();
                //hxy.print_bin_counts();
                
                double granularity_x = (double) (1.0/(double) x_num_bins); 
                double granularity_y = (double) (1.0/(double) y_num_bins); 
                
                double mi=0.0;
                double px, py, pxy, x, y;
        	for(x=0.0000001;x<1.0;x=x+granularity_x) {//Starting with x=0.0000001 to avoid double inaccuracy
                        px=hx.probability(x);
        		for(y=0.0000001;y<1.0;y=y+granularity_y) {
        			pxy=hxy.probability(x,y);
        			py=hy.probability(y);
                                if (px == 1 || py == 1) {//Diversity experiment (to ensure neuron value is not fixed)
                                        mi = 1;
                                }
                                else if(px!=0.0 && py!=0.0 && pxy !=0.0)
        				//mi+=pxy*((double)log2((double)(pxy/(px*py))));//*(0.05f*0.05f);
        				mi+=pxy*(log_anybase((double)(pxy/(px*py)), y_num_bins));//Assumption: x_num_bins = y_num_bins;
        		}
        	}
                if (mi > 1.0) {
                        std::cout<<" MUTUAL INFORMATION IS MORE THAN : 1.0 "<<mi<<std::endl;//Assuming x_num_bins = y_num_bins
                        exit(0);
                }
                //Scale mi between 0-1 by dividing by its maximum value
                //mi = mi/log2((double)y_num_bins); //Assuming x_num_bins = y_num_bins
        	return mi;		
        }

};

double compute_standev_scaled(int num_bin, int y_x_delay, int active_time_steps, std::vector <double> sequence_x) {

  //std::vector <double> X; //One list for each history step
  double standev, standev_max, standev_min, standev_scaled; 

  ////Join data from all the trials into a single vector list 
  //for (int i=0; i<sequence_x.size(); i++) {//For each trial
  //        for (int j=y_x_delay; j<=active_time_steps-1; j++) {
  //                X.push_back(sequence_x[i][j]);//Store X for each history step in a separate vector
  //        }
  //}

  standev = histogram::standard_deviation(sequence_x, num_bin);
  
  //Scale standard deviation to a value of 0-1
  standev_min = 0.0; 
  standev_max = (1.0 - 0.0)/2.0; // (X_max - X_min)/2, where X_max = 1 (Issue: Not true for ReLU) and X_min = 0
  standev_scaled = (standev - standev_min)/(standev_max - standev_min); //Scaling it to 0-1
  return standev_scaled;
}

double compute_entropy(int num_bin,std::vector <double> sequence_x) {

  double entropy; 

  entropy = histogram::entropy(sequence_x, num_bin);
  return entropy;
}

double compute_mutual_information(int num_bin, std::vector <double> sequence_x, std::vector <double> sequence_y) {
 
  double mutual_information = 0.0; 
  int x_num_bins = num_bin; 
  int y_num_bins = num_bin;

  mutual_information += histogram::mutual_inf(sequence_x, sequence_y, x_num_bins, y_num_bins); //
  return mutual_information;
}

void op_ip_mutual_info (std::vector< std::vector<double> > input_sequences, std::vector< std::vector< std::vector <double> > > output_sequences, int y_x_delay, int num_bin, int active_time_steps, int max_history, int min_history) {
          
          for (int i=0; i<=output_sequences.size()-1; i++) {//For each output node 
                  for (int j=0; j<=y_x_delay; j++) {//For each history step
                          std::vector<double> input, output;
                          for (int k = 0; k <= input_sequences.size()-1; k++) { //For each trial

                                  //Align input, output for different history steps from all the trials
                                  std::vector<double>::const_iterator ip_first = input_sequences[k].begin();
                                  std::vector<double>::const_iterator ip_last  = input_sequences[k].begin()+(input_sequences[k].size()) - j;
                                  std::vector<double> trial_input(ip_first, ip_last);
                                  input.insert(input.end(), trial_input.begin(), trial_input.end() );
                                  
                                  std::vector<double>::const_iterator op_first = output_sequences[i][k].begin() + j;
                                  std::vector<double>::const_iterator op_last  = output_sequences[i][k].end();
                                  std::vector<double> trial_output(op_first, op_last);
                                  output.insert(output.end(), trial_output.begin(), trial_output.end() );
                          }
                          if (input.size() != output.size()) {
                                  std::cout<<" ERROR: Input and Output size mismatch in op_ip_mutual_info()"<<std::endl;
                                  exit(0);
                          }

                          double temp = compute_mutual_information(num_bin, input, output);
                          std::cout<<"Mutual Information between Output node: "<<i<<" and Input step: "<<j<<" = "<<temp<<std::endl;
                          input.clear();
                          output.clear();
                  } 
          }
}

double normalized_mi(const vector <double> &X, const vector <double> &Y, const double entropy_X, const double entropy_Y, int K, int num_bin) {
        double large_number = 100.0;
        double mi, norm_mi;
        double max_entropy = 1.0; //This the fraser entropy for a completely random variable (size = 80x1, generated by MATLAB)
        if (entropy_X == 0.0 || entropy_Y == 0.0) {//If one of the outputs is fixed
                norm_mi = 1.0; //New Implementation:: To ensure that fitness > 0. Old Implementation: Large number (Avoid division by zero)
        }
        else {
                mi = compute_mutual_information(num_bin, X, Y);
                std::cout<<" Binned Entropy: "<<entropy_X<<" "<<entropy_Y<<std::endl;
                //mi = slow_kraskov_mutual_information(K, X, Y);
                //mi = fraser_mutual_information(X, Y);
                std::cout<<"Binned MI: "<<mi<<std::endl;
                if (entropy_X<entropy_Y) {
                        //if(entropy_X < 0.5) { //If entropy is very small, no point looking at the mutual information
                        //        norm_mi = (1.0/entropy_X);
                        //}
                        //else {
                                //norm_mi = mi/(entropy_X/max_entropy);  //ref: Pablo estez (2009) - Normalized mutual information
                                norm_mi = (mi + (max_entropy-entropy_X))/2;  //ref: Pablo estez (2009) - Normalized mutual information
                        //} 
                }
                else {
                        //if(entropy_Y < 0.5) {//If entropy is very small, no point looking at the mutual information
                        //        norm_mi = (1.0/entropy_Y);
                        //}
                        //else {
                                //norm_mi = mi/(entropy_Y/max_entropy);  //ref: Pablo estez (2009) - Normalized mutual information
                                norm_mi = (mi + (max_entropy-entropy_Y))/2;  //ref: Pablo estez (2009) - Normalized mutual information
                        //} 
                }
        }
        return norm_mi;
}


bool memory_evaluate(Organism *org, int generation, int org_index, int num_active_outputs, int output_start_index, int output_end_index, const std::vector < vector < double > > &independent_archive, const std::vector < vector <double> > &output_activations, int max_nw_size) {
  
  double errorsum = 0.0;
  int num_output_nodes = org->net->outputs.size();
  int num_input_nodes = org->net->inputs.size();
  double large_number = 100.0;
  int nw_size;

  //Parameters for information objective (Used to specify history window)
  int num_bin = 10; //Real value from 0-1 is discretized into these bins (Set to 2 for binary inputs)
  int K = 3; //Nearest neighbor parameter for Kraskov mutual information computation 
  //Print to file for plotting these parameters
  if (generation == 1) {
          std::cout<<" num_bin: "<<num_bin<<" Kraskov k: "<<K<<" num_outputs: "<<num_output_nodes<<std::endl;  
  }

  ////Write the features in output files
  //std::vector <double> row;
  //int j, k;
  //for ( j = 0 ; j <num_output_nodes ; j++ ) { 
  //        char temp[50];
  //        sprintf (temp, "output_features_%d", j+1+num_input_nodes);
  //        ofstream output_file(temp);
  //        ostream_iterator<double> output_iterator(output_file, " ");
  //        row.clear();
  //        for (k = 0; k < output_activations[0].size(); k++) {
  //                row.push_back(output_activations[j][k]);
  //                copy(row.begin(), row.end(), output_iterator);
  //                output_file  << '\n';
  //                row.clear(); 
  //        }
  //        output_file.close();
  //}

  //Fitness
  if (org->error !=1) { //If all outputs were activated

    double mutual_information = 0.0; //Mutual information range: 0 to Large number
    vector <double> entropy(output_end_index-output_start_index, 0.0);
    vector <double> entropy_archive(independent_archive.size(), 0.0);
    int mi_count = 0;
    int feature_count = 1;
    //Compute entropy of each output variable  
    for (int i=output_start_index; i<output_end_index; i++) {
             //entropy[i-output_start_index] = fraser_entropy(output_activations[org_index*num_active_outputs+i], feature_count); //Ranges between 0-1
             //feature_count++;
             entropy[i-output_start_index] = compute_entropy(num_bin, output_activations[org_index*num_active_outputs+i]); //Ranges between 0-1
             //std::cout<<"Entropy of Output Node ID: "<<org->net->outputs[independent_archive.size()+i]->node_id<<" : "<<entropy[i]<<std::endl;
    }     
    //Entropy of the Archived features (Future Work: No need to recompute this for the frozen network)
    for (int i=0; i<independent_archive.size(); i++) {
             //entropy_archive[i] = fraser_entropy(independent_archive[i], feature_count); //Ranges between 0-max_entropy
             //feature_count++;
             entropy_archive[i] = compute_entropy(num_bin, independent_archive[i]); //Ranges between 0-1
             //std::cout<<"Entropy of Archived Feature "<<i<<" : "<<entropy_archive[i]<<std::endl;
    }       

    //Compute pairwise NORMALIZED mutual information between all the newly added output nodes and between new outputs and the archive
    for (int i=output_start_index; i<output_end_index; i++) {
            double norm_mi;
            for (int j=output_start_index; j<output_end_index; j++) {
                    if (i != j && j > i) {//Skip duplicate pairs 
                             norm_mi = normalized_mi(output_activations[org_index*num_active_outputs+i], 
                                                             output_activations[org_index*num_active_outputs+j], entropy[i-output_start_index], entropy[j-output_start_index], K, 2);
                             if (norm_mi > mutual_information) { //Largest mutual information pair (No need to average)
                                             mutual_information = norm_mi;
                             }
                             //std::cout<<"Mutual Information between Outputs Node IDs: "<<org->net->outputs[independent_archive.size()+i]->node_id<<" "<<org->net->outputs[independent_archive.size()+j]->node_id<<": "<<norm_mi<<std::endl;
                    }
            }
            for (int j=0; j<independent_archive.size(); j++) {
                    norm_mi = normalized_mi(output_activations[org_index*num_active_outputs+i], 
                                                    independent_archive[j], entropy[i-output_start_index], entropy_archive[j], K, 2);
                    if (norm_mi > mutual_information) { //Largest mutual information pair (No need to average)
                                    mutual_information = norm_mi;
                    }
                    //std::cout<<"Mutual Information between Archived Feature and Output Node ID: "<<j<<" "<<org->net->outputs[independent_archive.size()+i]->node_id<<": "<<norm_mi<<std::endl;

            }
    }
    org->fitness1 = 0.0001 + (1-mutual_information)*100; //(((1-mutual_information) + entropy)/2)*100; //Minimize Mutual Info and Max variable entropy 
    //Fitness Penalty for Large networks
    nw_size = ((org->gnome)->compute_genome_size());
    //org->fitness2 = 100-genome_size; //(1.0/org->gnome->nodes.size())*100.0;//Cost for the size of the network
  }
  else {
    //The network is flawed (shouldnt happen)
    std::cout << " Net not activating. No path to output in genome id: "<<org->gnome->genome_id<<std::endl;//memory_startgenes makes sure no floating outputs
    errorsum=999.0;
    org->fitness1=0.0;//(1.0-large_number)*100;
    //org->fitness2=0.0;//Cost for the size of the network
  }

  if (org->fitness1>=99.000 && nw_size<=max_nw_size) {
        org->winner = true;
        std::cout<<"OBJECTIVE1 ACHIEVED:: WINNER FOUND"<<std::endl;
  }
  else{
        org->winner = false;
  }
  if (org->fitness2>=95.0) {
        //org->winner = true;
        std::cout<<"OBJECTIVE2 ACHIEVED:: WINNER FOUND"<<std::endl;
  }
  //if (org->fitness1>=99.0 && org->fitness2>=95.0) {
  //        org->winner = true;
  //}
  //else {
  //        org->winner = false;
  //}

  #ifndef NO_SCREEN_OUT
  cout<<"Org "<<(org->gnome)->genome_id<<"                                     error: "<<org->error<<endl;
  cout<<"Org "<<(org->gnome)->genome_id<<"                                     fitness1: "<<org->fitness1<<endl;
  cout<<"Org "<<(org->gnome)->genome_id<<"                                     fitness2: "<<org->fitness2<<endl;
  cout<<"Org "<<(org->gnome)->genome_id<<"                                     Size: "<<nw_size<<endl;
  #endif
  //exit(0);
  org->evaluated = true; //Aditya: for speed-up by preventing re-evaluation of the elites
  return org->winner;
}

//Print network input and outputs for each organism
void  print_ip_op(const std::vector < vector < double > > &input_data, const std::vector < vector <double> > &output_activations, int org_count, int num_active_outputs, int output_start_index, int output_end_index) {
        for (int org_index=0; org_index<org_count; org_index++) {
                std::cout<<"ORG NUMBER:: "<<org_index<<std::endl;
                for (int i=0; i<input_data.size(); i++) {
                        int output_count = 0; //For indexing the correct network output
                        //std::cout<<"Input: "<<input_data[i][0]<<" "<<input_data[i][1]<<" Output: ";  
                        for (int j=output_start_index; j<output_end_index; j++) {//Storing output from each non-frozen output nodes 
                                std::cout<<output_activations[org_index*num_active_outputs + output_count][i]<<"      ";
                                output_count++;
                        }
                        //std::cout<<std::endl;
                }
        } 
}

//Scale each output node values such that its value lie in the range of (0-0.9999)
void scale_output_activations(std::vector < vector <double> > &output_activations, int org_count, int num_active_outputs, int output_start_index, int output_end_index) {
        for (int org_index=0; org_index<org_count; org_index++) {
               int output_count = 0; //For indexing the correct network output  
               for (int i=output_start_index; i<output_end_index; i++) {//Storing output from each non-frozen output node
                      double max_value, min_value;
                      vector<double> single_output_activations = output_activations[org_index*num_active_outputs+output_count];
                      //min_value = *min_element(single_output_activations.begin(), single_output_activations.end()); 
                      //max_value = *max_element(single_output_activations.begin(), single_output_activations.end()); 
                      min_value = -1; //tanh activation for LSTM outputs
                      max_value = 1; //tanh activation for LSTM outputs
                      for (int j=0; j< single_output_activations.size(); j++) {
                              if (max_value==min_value) {
                                      single_output_activations[j] = 0.5; //Fix it to a constant
                              }
                              else {
                                      single_output_activations[j] = (single_output_activations[j] - min_value)/(max_value-min_value);
                              }
                              if (single_output_activations[j] == 1.0) {//Kraskov mutual information implementation cannot handle values >= 1.0
                                      single_output_activations[j] = single_output_activations[j]  - 0.0001;
                              }
                      }
                      output_activations[org_index*num_active_outputs+output_count] = single_output_activations;
                      output_count = output_count + 1;
               }
        }
}

void save_output_activations(int org_index, int num_active_outputs, std::vector < vector <double> > &output_activations, Network *net, int output_start_index, int output_end_index) {
        int output_count = 0; //For indexing the correct network output  
        for (int i=output_start_index; i<output_end_index; i++) {//Storing output from each non-frozen output node 
                //if ((net->outputs[i])->activation == 1.0){//Kraskov mutual information implementation cannot handles values >= 1.0
                //        (net->outputs[i])->activation = (net->outputs[i])->activation - 0.0000001;
                //}
                //std::cout<<((net->outputs[i])->activation)<<" ";;
                output_activations[org_index*num_active_outputs+output_count].push_back((net->outputs[i])->activation);
                output_count = output_count + 1;
        }
        //std::cout<<std::endl;
}

void memory_activate(Organism *org, int org_index, int num_active_outputs, int output_start_index, int output_end_index, const std::vector < vector < double > > &input_data, std::vector < vector <double> > &output_activations, int num_trials, int aisle_length) {
  
  Network *net;
  bool success;  //Check for successful activation
  int count=0;
  double output_error = 0.0;
  double output_error_1_step = 0.0;
  double output_error_2_step = 0.0;
  double average_output_error;
  double average_output_error_1_step;
  double average_output_error_2_step;
  double in[2]; //3-bit input - Bias, distance from T-junction, Binary light (older -> input data valid(1/0), Recall(1)/Instruct(0))
  in[0] = 1.0; //First input is Bias signal
  //int net_depth; //The max depth of the network to be activated
  //int num_trials = 20;
  //int max_rand_activate = 20;
  net=org->net;
  double reward = 0.0;
  std::vector <double> output_sequence;
  std::vector <double> output_sequence2;
  int trial_length = aisle_length*NEAT::input_sequence_len;//Total number of steps that an agent takes in a trial
  //net_depth=net->max_depth();
  
  //Load and activate the network on each input
  //for (int trial=0;trial< num_trials; trial++) {//Repeat several times for the same set of sequence
  for(int trial=0; trial < num_trials; trial++) {//For each new input sequence
           for (int step=0;step< trial_length; step++) { //WRITE PHASE:: For each step in the input sequence
  
                    //Activate NN
                    //in[1] = input_data[trial*trial_length+step][0]; //Input T-junction Distance
                    in[1] = input_data[trial*trial_length+step][1]; //Input Binary Light
                    
                    //in[3] = 0; //Recall signal
                    //std::cout<<"Input Sequence: ";
                    //std::cout<<in[1]<<":::";
                    net->load_sensors(in);
                    success=net->activate(); 
                    //std::cout<<std::endl<<" Some Output: "<<std::endl;
                    //std::cout<<" Some Output: "<<(net->outputs[0])->activation <<" "<<(net->outputs[1])->activation<<std::endl;
                    //save_output_activations(org_index, num_active_outputs, output_activations, net, output_start_index, output_end_index);
                    if (!success) {
                            org->error = 1;
                            std::cout << " Net not activating. No path to output in genome id: "<<org->gnome->genome_id<<std::endl;//memory_startgenes makes sure no floating outputs
                            break;//Breaks inner for-loop
                            //sprintf (temp, "not_activating_%d", org_index);
                            //std::ofstream outFile("not_activating",std::ios::out);
                            //((org)->gnome)->print_to_file(outFile);
                    }
                    save_output_activations(org_index, num_active_outputs, output_activations, net, output_start_index, output_end_index);
           }
           //if (org->error == 1){
           //         break;//Breaks outer for-loop
           //}
           //else { 
           //        for (int step=0;step< input_data[seqnum].size(); step++) {//RECALL PHASE:: For each step in the input sequence

           //                //Insert Zero (Don't-care) input activations at random time-steps
           //                in[1] = 0;//Input data
           //                in[2] = 0;//Input Data Valid
           //                in[3] = step+1;//Recall signal
           //                int rand_num = round(max_rand_activate*randfloat())+max_rand_activate;//ranges between 10-110
           //                //Activate NN for random time-steps
           //                //std::cout<<std::endl<<" Random Output: "<<std::endl;
           //                for (int i=0; i<rand_num; i++) {//Don't-care activate 
           //                        net->load_sensors(in); 
           //                        success=net->activate(); 
           //                        save_output_activations(org_index, num_active_outputs, output_activations, net, output_start_index, output_end_index);
           //                        //std::cout<<(net->outputs[0])->activation <<std::endl;
           //                }

           //                //Activate NN once after random-time steps for recall
           //                in[1] = 0;
           //                in[2] = 0;//Negative value indicates RECALL
           //                in[3] = step+1;//Recall signal
           //                net->load_sensors(in); //Give zeroes as input during recall phase
           //                success=net->activate();
           //                
           //                //std::cout<<std::endl<<" Actual Output: "<<std::endl;
           //                //std::cout<<(net->outputs[0])->activation <<std::endl;
           //                save_output_activations(org_index, num_active_outputs, output_activations, net, output_start_index, output_end_index);
           //                //***************START OLD DISCRETE INPUT****************
           //                ////First output is the actual stored value we are looking for 
           //                //if (((net->outputs[0])->activation >= 0.5 && input_data[seqnum][step]==-1.0) ||
           //                //    ((net->outputs[0])->activation < 0.5 && input_data[seqnum][step]==1.0)){
           //                //        output_error = output_error + 1;
           //                //}
           //                //count = count + 1;
           //                //***************END OLD DISCRETE INPUT****************

           //                //***************START NEW DISCRETE INPUT WITH SHAPING REWARDS****************
           //                //if ((net->outputs[0])->activation >= 0.5) {
           //                //        output_sequence.push_back(1.0);
           //                //}
           //                //else {
           //                //        output_sequence.push_back(-1.0);
           //                //}
           //                //***************END NEW DISCRETE INPUT WITH SHAPING REWARDS****************
           //                
           //                output_sequence.push_back((net->outputs[0])->activation);
           //                //output_sequence2.push_back((net->outputs[1])->activation);
           //                        
           //                output_error = output_error + abs((net->outputs[0])->activation - input_data[seqnum][step]);
           //                count = count + 1;

           //        }
           //        //if ((output_sequence[0]==input_data[seqnum][0]) && (output_sequence[1]==input_data[seqnum][1]) && (output_sequence[2]==input_data[seqnum][2])){
           //        //        reward = reward + 1;
           //        //}
           //        //else if ((output_sequence[0]==input_data[seqnum][0]) && (output_sequence[1]==input_data[seqnum][1])){
           //        //        reward = reward + 0.5;
           //        //}
           //        //if ((output_sequence[0]==input_data[seqnum][0])) {
           //        //        reward = reward + 1;
           //        //}
           //        //output_error_1_step = output_error_1_step + abs(output_sequence[0] - input_data[seqnum][2]);
           //        //output_error_2_step = output_error_2_step + abs(output_sequence2[0] - input_data[seqnum][1]);
           //        //if (abs(output_sequence[0] - input_data[seqnum][0]) < 0.05) {
           //        //        output_error = output_error + abs(output_sequence[1] - input_data[seqnum][1]);
           //        //}
           //        //else {
           //        //        output_error = output_error + abs(output_sequence[0] - input_data[seqnum][0]);
           //        //}
           //        //count = count+1;
           //        output_sequence.clear();
           //        output_sequence2.clear();
           //}
           net->flush(); //Flush after each sequence
  }
  //if (org->error == 1){
  //         break;//Breaks outer for-loop
  //}
  //}
  if (org->error == 1){
        ////The network is flawed (shouldnt happen)
        //std::cout << " Net not activating. No path to output in genome id: "<<org->gnome->genome_id<<std::endl;//memory_startgenes makes sure no floating outputs
        org->fitness1=0.0;
        //org->fitness2=0.0;
  }
  else {
        //average_output_error = output_error/count;//Ranges between 0-1
        //average_output_error_1_step = output_error_1_step/count;//Ranges between 0-1
        //average_output_error_2_step = output_error_2_step/count;//Ranges between 0-1
        //reward = 1.0-(average_output_error_1_step+average_output_error_2_step)/2; //1 Step recall reward 
        //reward = 1-average_output_error; //2 Step recall reward 

        //if (average_output_error_2_step < 0.05) {
        //        reward = 1-average_output_error_1_step; //2 Step recall reward 
        //}
        //else {
        //        reward = 0.5-average_output_error_2_step; //1 Step recall reward
        //}
        //org->fitness1 = (reward)*100;//(1-average_output_error)*100; //Ranges between 0-100  
        //org->fitness2 = org->fitness1; //std_dev*100; //To scale it to 0-100
  }
  //if (org->fitness1>=99.000) {
  //      //org->winner = true;
  //      std::cout<<"OBJECTIVE1 ACHIEVED:: WINNER FOUND"<<std::endl;
  //}
  //else{
  //     //org->winner = false;
  //}
  //if (org->fitness2>=99.0) {
  //      //org->winner = true;
  //      std::cout<<"OBJECTIVE2 ACHIEVED:: WINNER FOUND"<<std::endl;
  //}
  //if (org->fitness1>=99.0 && org->fitness2>=99.0) {
  //        org->winner = true;
  //}
  //else {
  //        org->winner = false;
  //}

  //#ifndef NO_SCREEN_OUT
  //cout<<"Org "<<(org->gnome)->genome_id<<"                                     error: "<<org->error<<endl;
  //cout<<"Org "<<(org->gnome)->genome_id<<"                                     fitness1: "<<org->fitness1<<endl;
  //cout<<"Org "<<(org->gnome)->genome_id<<"                                     fitness2: "<<org->fitness2<<endl;
  //#endif
  ////exit(0);
  //org->evaluated = true; //Aditya: for speed-up by preventing re-evaluation of the elites
  //return org->winner;

}

int memory_epoch(Population *pop,int generation,char *filename,int &winnernum,int &winnergenes,int &winnernodes, const std::vector < vector < double > > &input_data, std::vector< vector <double> > &independent_archive, int num_trials, int aisle_length, int max_nw_size) {
  vector<Organism*>::iterator curorg;
  vector<Species*>::iterator curspecies;
  bool win=false;
  int org_count = 0;
  int num_output_nodes = pop->organisms[0]->net->outputs.size();
  int output_start_index;
  int output_end_index;
  int num_active_outputs; 
  std::vector<int> unevaluated_org; //Stores the organism number for the individuals which haven't been evaluated in the population 
  time_t start_seconds, end_seconds;
  int genome_size;
  int new_hidden_nodes; 

  //If the startgenome was frozen, create the archive of the independent activations (of frozen outputs) at the start
  if (NEAT::frozen_startgenome == 1 && generation==1) {
          num_active_outputs = num_output_nodes-NEAT::batch_size;
          output_start_index = 0;
          output_end_index = num_active_outputs;
          independent_archive.resize(num_active_outputs, std::vector<double>(input_data.size()));
          //std::cout<<" num_active_outputs: "<< num_active_outputs<<std::endl;
          for (int i=0; i < pop->organisms.size(); i++) {
                  memory_activate(pop->organisms[i], 0, num_active_outputs, output_start_index, output_end_index, input_data, independent_archive, num_trials, aisle_length);
                  if (pop->organisms[i]->error !=1) { //If all outputs were activated
                          break;
                  }
          }
          //Scale output values to a range (0-0.999)
          scale_output_activations(independent_archive, 1, num_active_outputs, output_start_index, output_end_index);
  }

  if (independent_archive.size()==0) { //If we are starting from scratch, activate all the output neurons
          num_active_outputs = num_output_nodes;
  }
  else {                               //If we already have stored some features (frozen), then consider only the non-frozen outputs 
          num_active_outputs = NEAT::batch_size;
  }
  
  //Store the organism number for the individuals which haven't been evaluated in the population 
  for (int i=0; i < pop->organisms.size(); i++) {
          if (pop->organisms[i]->evaluated==false) {
                  unevaluated_org.push_back(i); 
          }
  }
  //int num_trials = 20;
  std::vector< std::vector <double> > output_activations(unevaluated_org.size()*num_active_outputs);//Vector of non-frozen output activations for each organism
  std::vector<int> org_win(unevaluated_org.size(), 0);//Win/Loss status for an organism

  //Activate networks of each non-evaluated organism with input data 
  output_start_index = num_output_nodes-num_active_outputs;
  output_end_index = num_output_nodes;
  #pragma omp parallel for //Parallelization of for loop 
  for (int i=0; i < unevaluated_org.size(); i++) {
          memory_activate(pop->organisms[unevaluated_org[i]], i, num_active_outputs, output_start_index, output_end_index, input_data, output_activations, num_trials, aisle_length);
  }
  
  //print_ip_op(input_data, output_activations, unevaluated_org.size(), num_active_outputs, output_start_index, output_end_index);
  //Scale output values to a range (0-0.999)
  scale_output_activations(output_activations, unevaluated_org.size(), num_active_outputs, output_start_index, output_end_index);

  //Print the network input and output for each organism
  //print_ip_op(input_data, output_activations, unevaluated_org.size(), num_active_outputs, output_start_index, output_end_index);
  
  //Evaluate each organism for its independent output features (number of features = num_active_outputs)
  output_start_index = 0; //Older version was --> First output value is used only for comparison with input. Remaining outputs are checked for independence
  output_end_index = num_active_outputs;
  //start_seconds = time(NULL); //Current Time in seconds since January 1, 1970
  #pragma omp parallel for //Parallelization of for loop 
  for (int i=0; i < unevaluated_org.size(); i++) {
          org_win[i] = memory_evaluate(pop->organisms[unevaluated_org[i]], generation, i, num_active_outputs, output_start_index, output_end_index, independent_archive, output_activations, max_nw_size);
  }
  //end_seconds = time(NULL); //Current Time in seconds since January 1, 1970
  //std::cout << "Total Organism Evaluation Time: "<< (double)(end_seconds-start_seconds)<< " seconds." << "\n";
 
  //Save the first winner 
  for (int i=0; i < unevaluated_org.size(); i++) {
          if (org_win[i]==1) {
            win=true;
            winnernum=pop->organisms[unevaluated_org[i]]->gnome->genome_id;
            winnergenes=pop->organisms[unevaluated_org[i]]->gnome->extrons();
            winnernodes=(pop->organisms[unevaluated_org[i]]->gnome->nodes).size();
            ////Push in the num_active_outputs features
            //for (int j=0; j<num_active_outputs; j++) {
            //        independent_archive.push_back(output_activations[i*num_active_outputs+j]);
            //}
            break;
          }
  }
  
  ////Average and max their fitnesses for dumping to file and snapshot
  //for(curspecies=(pop->species).begin();curspecies!=(pop->species).end();++curspecies) {

  //  //This experiment control routine issues commands to collect ave
  //  //and max fitness, as opposed to having the snapshot do it, 
  //  //because this allows flexibility in terms of what time
  //  //to observe fitnesses at

  //  (*curspecies)->compute_average_fitness();
  //  (*curspecies)->compute_max_fitness();
  //}

  //Take a snapshot of the population, so that it can be
  //visualized later on
  //if ((generation%1)==0)
  //  pop->snapshot();

  //Only print to file every print_every generations
  if  (win) {//||
//       ((generation%(NEAT::print_every))==0)) //Print every generation happens inside epoch_multiobj
    //pop->print_to_file_by_species(filename);
  }

  if (win) {
    for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
      if ((*curorg)->winner) {
        genome_size = ((*curorg)->gnome)->compute_genome_size();
        new_hidden_nodes = ((*curorg)->gnome)->nodes.size()-2-2; //When starting with two input nodes and two output nodes 
	cout<<"WINNER IS #"<<((*curorg)->gnome)->genome_id<<" Added Hidden Nodes "<<new_hidden_nodes <<" Generation IS "<<generation<<" Num outputs "<<((*curorg)->net->outputs.size()) <<" Size IS "<< genome_size<<endl;
	//Prints the winner to file
	//IMPORTANT: This causes generational file output!
	char temp[200];
	sprintf (temp, "/scratch/cluster/aditya/memory_expt_files/winner/winner_%d_%d", generation,((*curorg)->net->outputs.size()));
	print_Genome_tofile((*curorg)->gnome,temp);
        break;
      }
    }
  }
  else {
//  if(generation <= 999)
        pop->epoch(generation, filename);
        //pop->epoch_multiobj(generation, filename); //Aditya (NSGA-2)
  //if  (win||
  //     ((generation%(NEAT::print_every))==0))
  //  pop->print_to_file_by_species(strcat(filename, "_del"));
  }

  if (win) return 1;
  else return 0;

}


//Perform evolution on XOR, for gens generations
Population *xor_test(int gens) {
    Population *pop=0;
    Genome *start_genome;
    char curword[20];
    int id;

    ostringstream *fnamebuf;
    int gen;
 
    int evals[NEAT::num_runs];  //Hold records for each run
    int genes[NEAT::num_runs];
    int nodes[NEAT::num_runs];
    int winnernum;
    int winnergenes;
    int winnernodes;
    //For averaging
    int totalevals=0;
    int totalgenes=0;
    int totalnodes=0;
    int expcount;
    int samples;  //For averaging

    memset (evals, 0, NEAT::num_runs * sizeof(int));
    memset (genes, 0, NEAT::num_runs * sizeof(int));
    memset (nodes, 0, NEAT::num_runs * sizeof(int));

    ifstream iFile("xorstartgenes",ios::in);

    cout<<"START XOR TEST"<<endl;

    cout<<"Reading in the start genome"<<endl;
    //Read in the start Genome
    iFile>>curword;
    iFile>>id;
    cout<<"Reading in Genome id "<<id<<endl;
    start_genome=new Genome(id,iFile);
    iFile.close();

    for(expcount=0;expcount<NEAT::num_runs;expcount++) {
      //Spawn the Population
      cout<<"Spawning Population off Genome2"<<endl;

      pop=new Population(start_genome,NEAT::pop_size);
      
      cout<<"Verifying Spawned Pop"<<endl;
      pop->verify();
      
      for (gen=1;gen<=gens;gen++) {
	cout<<"Epoch "<<gen<<endl;	

	//This is how to make a custom filename
	fnamebuf=new ostringstream();
	(*fnamebuf)<<"gen_"<<gen<<ends;  //needs end marker

	#ifndef NO_SCREEN_OUT
	cout<<"name of fname: "<<fnamebuf->str()<<endl;
	#endif

	char temp[50];
	sprintf (temp, "gen_%d", gen);

	//Check for success
	if (xor_epoch(pop,gen,temp,winnernum,winnergenes,winnernodes)) {
	  //	if (xor_epoch(pop,gen,fnamebuf->str(),winnernum,winnergenes,winnernodes)) {
	  //Collect Stats on end of experiment
	  evals[expcount]=NEAT::pop_size*(gen-1)+winnernum;
	  genes[expcount]=winnergenes;
	  nodes[expcount]=winnernodes;
	  gen=gens;

	}
	
	//Clear output filename
	fnamebuf->clear();
	delete fnamebuf;
	
      }

      if (expcount<NEAT::num_runs-1) delete pop;
      
    }

    //Average and print stats
    cout<<"Nodes: "<<endl;
    for(expcount=0;expcount<NEAT::num_runs;expcount++) {
      cout<<nodes[expcount]<<endl;
      totalnodes+=nodes[expcount];
    }
    
    cout<<"Genes: "<<endl;
    for(expcount=0;expcount<NEAT::num_runs;expcount++) {
      cout<<genes[expcount]<<endl;
      totalgenes+=genes[expcount];
    }
    
    cout<<"Evals "<<endl;
    samples=0;
    for(expcount=0;expcount<NEAT::num_runs;expcount++) {
      cout<<evals[expcount]<<endl;
      if (evals[expcount]>0)
      {
	totalevals+=evals[expcount];
	samples++;
      }
    }

    cout<<"Failures: "<<(NEAT::num_runs-samples)<<" out of "<<NEAT::num_runs<<" runs"<<endl;
    cout<<"Average Nodes: "<<(samples>0 ? (double)totalnodes/samples : 0)<<endl;
    cout<<"Average Genes: "<<(samples>0 ? (double)totalgenes/samples : 0)<<endl;
    cout<<"Average Evals: "<<(samples>0 ? (double)totalevals/samples : 0)<<endl;

    return pop;

}

bool xor_evaluate(Organism *org) {
  Network *net;
  double out[4]; //The four outputs
  double this_out; //The current output
  int count;
  double errorsum;

  bool success;  //Check for successful activation
  int numnodes;  /* Used to figure out how many nodes
		    should be visited during activation */

  int net_depth; //The max depth of the network to be activated
  int relax; //Activates until relaxation

  //The four possible input combinations to xor
  //The first number is for biasing
  double in[4][3]={{1.0,0.0,0.0},
		   {1.0,0.0,1.0},
		   {1.0,1.0,0.0},
		   {1.0,1.0,1.0}};
  
  net=org->net;
  numnodes=((org->gnome)->nodes).size();

  net_depth=net->max_depth();

  //TEST CODE: REMOVE
  //cout<<"ACTIVATING: "<<org->gnome<<endl;
  //cout<<"DEPTH: "<<net_depth<<endl;

  //Load and activate the network on each input
  for(count=0;count<=3;count++) {
    net->load_sensors(in[count]);

    //Relax net and get output
    success=net->activate();

    //use depth to ensure relaxation
    for (relax=0;relax<=net_depth;relax++) {
      success=net->activate();
      this_out=(*(net->outputs.begin()))->activation;
    }

    out[count]=(*(net->outputs.begin()))->activation;

    net->flush();

  }
  
  if (success) {
    errorsum=(fabs(out[0])+fabs(1.0-out[1])+fabs(1.0-out[2])+fabs(out[3]));
    org->fitness1=pow((4.0-errorsum),2);
    org->error=errorsum;
  }
  else {
    //The network is flawed (shouldnt happen)
    errorsum=999.0;
    org->fitness1=0.001;
  }

  #ifndef NO_SCREEN_OUT
  cout<<"Org "<<(org->gnome)->genome_id<<"                                     error: "<<errorsum<<"  ["<<out[0]<<" "<<out[1]<<" "<<out[2]<<" "<<out[3]<<"]"<<endl;
  cout<<"Org "<<(org->gnome)->genome_id<<"                                     fitness: "<<org->fitness1<<endl;
  #endif

  //  if (errorsum<0.05) { 
  //if (errorsum<0.2) {
  if ((out[0]<0.5)&&(out[1]>=0.5)&&(out[2]>=0.5)&&(out[3]<0.5)) {
    org->winner=true;
    return true;
  }
  else {
    org->winner=false;
    return false;
  }

}

int xor_epoch(Population *pop,int generation,char *filename,int &winnernum,int &winnergenes,int &winnernodes) {
  vector<Organism*>::iterator curorg;
  vector<Species*>::iterator curspecies;
  //char cfilename[100];
  //strncpy( cfilename, filename.c_str(), 100 );

  //ofstream cfilename(filename.c_str());

  bool win=false;


  //Evaluate each organism on a test
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
    if (xor_evaluate(*curorg)) {
      win=true;
      winnernum=(*curorg)->gnome->genome_id;
      winnergenes=(*curorg)->gnome->extrons();
      winnernodes=((*curorg)->gnome->nodes).size();
      if (winnernodes==5) {
	//You could dump out optimal genomes here if desired
	//(*curorg)->gnome->print_to_filename("xor_optimal");
	//cout<<"DUMPED OPTIMAL"<<endl;
      }
    }
  }
  
  //Average and max their fitnesses for dumping to file and snapshot
  for(curspecies=(pop->species).begin();curspecies!=(pop->species).end();++curspecies) {

    //This experiment control routine issues commands to collect ave
    //and max fitness, as opposed to having the snapshot do it, 
    //because this allows flexibility in terms of what time
    //to observe fitnesses at

    (*curspecies)->compute_average_fitness();
    (*curspecies)->compute_max_fitness();
  }

  //Take a snapshot of the population, so that it can be
  //visualized later on
  //if ((generation%1)==0)
  //  pop->snapshot();

  //Only print to file every print_every generations
  if  (win||
       ((generation%(NEAT::print_every))==0))
    pop->print_to_file_by_species(filename);


  if (win) {
    for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
      if ((*curorg)->winner) {
	cout<<"WINNER IS #"<<((*curorg)->gnome)->genome_id<<endl;
	//Prints the winner to file
	//IMPORTANT: This causes generational file output!
	print_Genome_tofile((*curorg)->gnome,"xor_winner");
      }
    }
    
  }

  pop->epoch(generation, filename);

  if (win) return 1;
  else return 0;

}

//Perform evolution on single pole balacing, for gens generations
Population *pole1_test(int gens) {
    Population *pop=0;
    Genome *start_genome;
    char curword[20];
    int id;

    ostringstream *fnamebuf;
    int gen;

    int expcount;
    int status;
    int runs[NEAT::num_runs];
    int totalevals;
    int samples;  //For averaging

    memset (runs, 0, NEAT::num_runs * sizeof(int));

    ifstream iFile("pole1startgenes",ios::in);

    cout<<"START SINGLE POLE BALANCING EVOLUTION"<<endl;

    cout<<"Reading in the start genome"<<endl;
    //Read in the start Genome
    iFile>>curword;
    iFile>>id;
    cout<<"Reading in Genome id "<<id<<endl;
    start_genome=new Genome(id,iFile);
    iFile.close();
  
    //Run multiple experiments
    for(expcount=0;expcount<NEAT::num_runs;expcount++) {

      cout<<"EXPERIMENT #"<<expcount<<endl;

      cout<<"Start Genome: "<<start_genome<<endl;
      
      //Spawn the Population
      cout<<"Spawning Population off Genome"<<endl;
      
      pop=new Population(start_genome,NEAT::pop_size);
      
      cout<<"Verifying Spawned Pop"<<endl;
      pop->verify();

      for (gen=1;gen<=gens;gen++) {
	cout<<"Generation "<<gen<<endl;
	
	fnamebuf=new ostringstream();
	(*fnamebuf)<<"gen_"<<gen<<ends;  //needs end marker

#ifndef NO_SCREEN_OUT
	cout<<"name of fname: "<<fnamebuf->str()<<endl;
#endif	

	char temp[50];
        sprintf (temp, "gen_%d", gen);

	status=pole1_epoch(pop,gen,temp);
	//status=(pole1_epoch(pop,gen,fnamebuf->str()));
	
	if (status) {
	  runs[expcount]=status;
	  gen=gens+1;
	}
	
	fnamebuf->clear();
	delete fnamebuf;
	
      }

      if (expcount<NEAT::num_runs-1) delete pop;
    }

    totalevals=0;
    samples=0;
    for(expcount=0;expcount<NEAT::num_runs;expcount++) {
      cout<<runs[expcount]<<endl;
      if (runs[expcount]>0)
      {
        totalevals+=runs[expcount];
        samples++;
      }
    }

    cout<<"Failures: "<<(NEAT::num_runs-samples)<<" out of "<<NEAT::num_runs<<" runs"<<endl;
    cout<<"Average evals: "<<(samples>0 ? (double)totalevals/samples : 0)<<endl;

    return pop;

}

int pole1_epoch(Population *pop,int generation,char *filename) {
  vector<Organism*>::iterator curorg;
  vector<Species*>::iterator curspecies;
  //char cfilename[100];
  //strncpy( cfilename, filename.c_str(), 100 );

  //ofstream cfilename(filename.c_str());

  bool win=false;
  int winnernum;

  //Evaluate each organism on a test
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
    if (pole1_evaluate(*curorg)) win=true;
  }
  
  //Average and max their fitnesses for dumping to file and snapshot
  for(curspecies=(pop->species).begin();curspecies!=(pop->species).end();++curspecies) {

    //This experiment control routine issues commands to collect ave
    //and max fitness, as opposed to having the snapshot do it, 
    //because this allows flexibility in terms of what time
    //to observe fitnesses at

    (*curspecies)->compute_average_fitness();
    (*curspecies)->compute_max_fitness();
  }

  //Take a snapshot of the population, so that it can be
  //visualized later on
  //if ((generation%1)==0)
  //  pop->snapshot();

  //Only print to file every print_every generations
  if  (win||
       ((generation%(NEAT::print_every))==0))
    pop->print_to_file_by_species(filename);

  if (win) {
    for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
      if ((*curorg)->winner) {
	winnernum=((*curorg)->gnome)->genome_id;
	cout<<"WINNER IS #"<<((*curorg)->gnome)->genome_id<<endl;
      }
    }    
  }

  //Create the next generation
  pop->epoch(generation, filename);

  if (win) return ((generation-1)*NEAT::pop_size+winnernum);
  else return 0;

}

bool pole1_evaluate(Organism *org) {
  Network *net;

  int numnodes;  /* Used to figure out how many nodes
		    should be visited during activation */
  int thresh;  /* How many visits will be allowed before giving up 
		  (for loop detection) */

  //  int MAX_STEPS=120000;
 int MAX_STEPS=100000;
  
  net=org->net;
  numnodes=((org->gnome)->nodes).size();
  thresh=numnodes*2;  //Max number of visits allowed per activation
  
  //Try to balance a pole now
  org->fitness1 = go_cart(net,MAX_STEPS,thresh);

#ifndef NO_SCREEN_OUT
  cout<<"Org "<<(org->gnome)->genome_id<<" fitness: "<<org->fitness1<<endl;
#endif

  //Decide if its a winner
  if (org->fitness1>=MAX_STEPS) { 
    org->winner=true;
    return true;
  }
  else {
    org->winner=false;
    return false;
  }

}

//     cart_and_pole() was take directly from the pole simulator written
//     by Richard Sutton and Charles Anderson.
int go_cart(Network *net,int max_steps,int thresh)
{
   float x,			/* cart position, meters */
         x_dot,			/* cart velocity */
         theta,			/* pole angle, radians */
         theta_dot;		/* pole angular velocity */
   int steps=0,y;

   int random_start=1;

   double in[5];  //Input loading array

   double out1;
   double out2;

//     double one_degree= 0.0174532;	/* 2pi/360 */
//     double six_degrees=0.1047192;
   double twelve_degrees=0.2094384;
//     double thirty_six_degrees= 0.628329;
//     double fifty_degrees=0.87266;

   vector<NNode*>::iterator out_iter;

   if (random_start) {
     /*set up random start state*/
     x = (lrand48()%4800)/1000.0 - 2.4;
     x_dot = (lrand48()%2000)/1000.0 - 1;
     theta = (lrand48()%400)/1000.0 - .2;
     theta_dot = (lrand48()%3000)/1000.0 - 1.5;
    }
   else 
     x = x_dot = theta = theta_dot = 0.0;
   
   /*--- Iterate through the action-learn loop. ---*/
   while (steps++ < max_steps)
     {
       
       /*-- setup the input layer based on the four iputs --*/
       //setup_input(net,x,x_dot,theta,theta_dot);
       in[0]=1.0;  //Bias
       in[1]=(x + 2.4) / 4.8;;
       in[2]=(x_dot + .75) / 1.5;
       in[3]=(theta + twelve_degrees) / .41;
       in[4]=(theta_dot + 1.0) / 2.0;
       net->load_sensors(in);

       //activate_net(net);   /*-- activate the network based on the input --*/
       //Activate the net
       //If it loops, exit returning only fitness of 1 step
       if (!(net->activate())) return 1;

      /*-- decide which way to push via which output unit is greater --*/
       out_iter=net->outputs.begin();
       out1=(*out_iter)->activation;
       ++out_iter;
       out2=(*out_iter)->activation;
       if (out1 > out2)
	 y = 0;
       else
	 y = 1;
       
       /*--- Apply action to the simulated cart-pole ---*/
       cart_pole(y, &x, &x_dot, &theta, &theta_dot);
       
       /*--- Check for failure.  If so, return steps ---*/
       if (x < -2.4 || x > 2.4  || theta < -twelve_degrees ||
	   theta > twelve_degrees) 
         return steps;             
     }
   
   return steps;
} 


//     cart_and_pole() was take directly from the pole simulator written
//     by Richard Sutton and Charles Anderson.
//     This simulator uses normalized, continous inputs instead of 
//    discretizing the input space.
/*----------------------------------------------------------------------
   cart_pole:  Takes an action (0 or 1) and the current values of the
 four state variables and updates their values by estimating the state
 TAU seconds later.
----------------------------------------------------------------------*/
void cart_pole(int action, float *x,float *x_dot, float *theta, float *theta_dot) {
  float xacc,thetaacc,force,costheta,sintheta,temp;
  
  const float GRAVITY=9.8;
  const float MASSCART=1.0;
  const float MASSPOLE=0.1;
  const float TOTAL_MASS=(MASSPOLE + MASSCART);
  const float LENGTH=0.5;	  /* actually half the pole's length */
  const float POLEMASS_LENGTH=(MASSPOLE * LENGTH);
  const float FORCE_MAG=10.0;
  const float TAU=0.02;	  /* seconds between state updates */
  const float FOURTHIRDS=1.3333333333333;

  force = (action>0)? FORCE_MAG : -FORCE_MAG;
  costheta = cos(*theta);
  sintheta = sin(*theta);
  
  temp = (force + POLEMASS_LENGTH * *theta_dot * *theta_dot * sintheta)
    / TOTAL_MASS;
  
  thetaacc = (GRAVITY * sintheta - costheta* temp)
    / (LENGTH * (FOURTHIRDS - MASSPOLE * costheta * costheta
		 / TOTAL_MASS));
  
  xacc  = temp - POLEMASS_LENGTH * thetaacc* costheta / TOTAL_MASS;
  
  /*** Update the four state variables, using Euler's method. ***/
  
  *x  += TAU * *x_dot;
  *x_dot += TAU * xacc;
  *theta += TAU * *theta_dot;
  *theta_dot += TAU * thetaacc;
}

/* ------------------------------------------------------------------ */
/* Double pole balacing                                               */
/* ------------------------------------------------------------------ */

//Perform evolution on double pole balacing, for gens generations
//If velocity is false, then velocity information will be withheld from the 
//network population (non-Markov)
Population *pole2_test(int gens,int velocity) {
    Population *pop=0;
    Genome *start_genome;
    char curword[20];
    int id;

    ostringstream *fnamebuf;
    int gen;
    CartPole *thecart;

    //Stat collection variables
    int highscore;
    int record[NEAT::num_runs][1000];
    double recordave[1000];
    int genesrec[NEAT::num_runs][1000];
    double genesave[1000];
    int nodesrec[NEAT::num_runs][1000];
    double nodesave[1000];
    int winnergens[NEAT::num_runs];
    int initcount;
    int champg, champn, winnernum;  //Record number of genes and nodes in champ
    int run;
    int curtotal; //For averaging
    int samples;  //For averaging

    ofstream oFile("statout",ios::out);

    champg=0;
    champn=0;

    //Initialize the stat recording arrays
    for (initcount=0;initcount<gens;initcount++) {
      recordave[initcount]=0;
      genesave[initcount]=0;
      nodesave[initcount]=0;
      for (run=0;run<NEAT::num_runs;++run) {
	record[run][initcount]=0;
	genesrec[run][initcount]=0;
	nodesrec[run][initcount]=0;
      }
    }
    memset (winnergens, 0, NEAT::num_runs * sizeof(int));

    char *non_markov_starter="pole2startgenes2";
    char *markov_starter="pole2startgenes1";
    char *startstring;

    if (velocity==0) startstring=non_markov_starter;
    else if (velocity==1) startstring=markov_starter;
    ifstream iFile(startstring,ios::in);
    //ifstream iFile("pole2startgenes",ios::in);

    cout<<"START DOUBLE POLE BALANCING EVOLUTION"<<endl;
    if (!velocity)
      cout<<"NO VELOCITY INPUT"<<endl;

    cout<<"Reading in the start genome"<<endl;
    //Read in the start Genome
    iFile>>curword;
    iFile>>id;
    cout<<"Reading in Genome id "<<id<<endl;
    start_genome=new Genome(id,iFile);
    iFile.close();

    cout<<"Start Genome: "<<start_genome<<endl;

    for (run=0;run<NEAT::num_runs;run++) {
      
      cout<<"RUN #"<<run<<endl;

      //Spawn the Population from starter gene
      cout<<"Spawning Population off Genome"<<endl;
      pop=new Population(start_genome,NEAT::pop_size);
      
      //Alternative way to start off of randomly connected genomes
      //pop=new Population(pop_size,7,1,10,false,0.3);

      cout<<"Verifying Spawned Pop"<<endl;
      pop->verify();
      
      //Create the Cart
      thecart=new CartPole(true,velocity);
      
      for (gen=1;gen<=gens;gen++) {
	cout<<"Epoch "<<gen<<endl;
	
	fnamebuf=new ostringstream();
	(*fnamebuf)<<"gen_"<<gen<<ends;  //needs end marker
#ifndef NO_SCREEN_OUT
	cout<<"name of fname: "<<fnamebuf->str()<<endl;
#endif

	char temp[50];
        sprintf (temp, "gen_%d", gen);

	highscore=pole2_epoch(pop,gen,temp,velocity, thecart,champg,champn,winnernum,oFile);
	//highscore=pole2_epoch(pop,gen,fnamebuf->str(),velocity, thecart,champg,champn,winnernum,oFile);  
	
	//cout<<"GOT HIGHSCORE FOR GEN "<<gen<<": "<<highscore-1<<endl;
	
	record[run][gen-1]=highscore-1;
	genesrec[run][gen-1]=champg;
	nodesrec[run][gen-1]=champn;
	
	fnamebuf->clear();
	delete fnamebuf;
	
	//Stop right at the winnergen
	if (((pop->winnergen)!=0)&&(gen==(pop->winnergen))) {
	  winnergens[run]=NEAT::pop_size*(gen-1)+winnernum;
	  gen=gens+1;
	}
	
	//In non-MARKOV, stop right at winning (could go beyond if desired)
	if ((!(thecart->MARKOV))&&((pop->winnergen)!=0))
	  gen=gens+1;

#ifndef NO_SCREEN_OUT
      cout<<"gen = "<<gen<<" gens = "<<gens<<endl;
#endif

      if (gen==(gens-1)) oFile<<"FAIL: Last gen on run "<<run<<endl;
      

      }

      if (run<NEAT::num_runs-1) delete pop;
      delete thecart;

    }

    cout<<"Generation highs: "<<endl;
    oFile<<"Generation highs: "<<endl;
    for(gen=0;gen<=gens-1;gen++) {
      curtotal=0;
      for (run=0;run<NEAT::num_runs;++run) {
	if (record[run][gen]>0) {
	  cout<<setw(8)<<record[run][gen]<<" ";
	  oFile<<setw(8)<<record[run][gen]<<" ";
	  curtotal+=record[run][gen];
	}
	else {
	  cout<<"         ";
	  oFile<<"         ";
	  curtotal+=100000;
	}
	recordave[gen]=(double) curtotal/NEAT::num_runs;
	
      }
      cout<<endl;
      oFile<<endl;
    }

    cout<<"Generation genes in champ: "<<endl;
    for(gen=0;gen<=gens-1;gen++) {
      curtotal=0;
      samples=0;
      for (run=0;run<NEAT::num_runs;++run) {
	if (genesrec[run][gen]>0) {
	  cout<<setw(4)<<genesrec[run][gen]<<" ";
	  oFile<<setw(4)<<genesrec[run][gen]<<" ";
	  curtotal+=genesrec[run][gen];
	  samples++;
	}
	else {
	  cout<<setw(4)<<"     ";
	  oFile<<setw(4)<<"     ";
	}
      }
      genesave[gen]=(double) curtotal/samples;

      cout<<endl;
      oFile<<endl;
    }

    cout<<"Generation nodes in champ: "<<endl;
    oFile<<"Generation nodes in champ: "<<endl;
    for(gen=0;gen<=gens-1;gen++) {
      curtotal=0;
      samples=0;
      for (run=0;run<NEAT::num_runs;++run) {
	if (nodesrec[run][gen]>0) {
	  cout<<setw(4)<<nodesrec[run][gen]<<" ";
	  oFile<<setw(4)<<nodesrec[run][gen]<<" ";
	  curtotal+=nodesrec[run][gen];
	  samples++;
	}
	else {
	  cout<<setw(4)<<"     ";
	  oFile<<setw(4)<<"     ";
	}
      }
      nodesave[gen]=(double) curtotal/samples;

      cout<<endl;
      oFile<<endl;
    }

    cout<<"Generational record fitness averages: "<<endl;
    oFile<<"Generational record fitness averages: "<<endl;
    for(gen=0;gen<gens-1;gen++) {
      cout<<recordave[gen]<<endl;
      oFile<<recordave[gen]<<endl;
    }

    cout<<"Generational number of genes in champ averages: "<<endl;
    oFile<<"Generational number of genes in champ averages: "<<endl;
    for(gen=0;gen<gens-1;gen++) {
      cout<<genesave[gen]<<endl;
      oFile<<genesave[gen]<<endl;
    }

    cout<<"Generational number of nodes in champ averages: "<<endl;
    oFile<<"Generational number of nodes in champ averages: "<<endl;
    for(gen=0;gen<gens-1;gen++) {
      cout<<nodesave[gen]<<endl;
      oFile<<nodesave[gen]<<endl;
    }

    cout<<"Winner evals: "<<endl;
    oFile<<"Winner evals: "<<endl;
    curtotal=0;
    samples=0;
    for (run=0;run<NEAT::num_runs;++run) {
      cout<<winnergens[run]<<endl;
      oFile<<winnergens[run]<<endl;
      if (winnergens[run]>0)
      {
        curtotal+=winnergens[run];
        samples++;
      }
    }
    cout<<"Failures: "<<(NEAT::num_runs-samples)<<" out of "<<NEAT::num_runs<<" runs"<<endl;
    oFile<<"Failures: "<<(NEAT::num_runs-samples)<<" out of "<<NEAT::num_runs<<" runs"<<endl;

    cout<<"Average # evals: "<<(samples>0 ? (double) curtotal/samples : 0)<<endl;
    oFile<<"Average # evals: "<<(samples>0 ? (double) curtotal/samples : 0)<<endl;

    oFile.close();

    return pop;

}

//This is used for list sorting of Species by fitness of best organism
//highest fitness first
//Used to choose which organism to test
//bool order_new_species(Species *x, Species *y) {
//
//  return (x->compute_max_fitness() > 
//	  y->compute_max_fitness());
//}

int pole2_epoch(Population *pop,int generation,char *filename,bool velocity,
		CartPole *thecart,int &champgenes,int &champnodes,
		int &winnernum, ofstream &oFile) {
  //char cfilename[100];
  //strncpy( cfilename, filename.c_str(), 100 );

  //ofstream cfilename(filename.c_str());

  vector<Organism*>::iterator curorg;
  vector<Species*>::iterator curspecies;

  vector<Species*> sorted_species;  //Species sorted by max fit org in Species

  int pause;
  bool win=false;

  double champ_fitness;
  Organism *champ;

  //double statevals[5]={-0.9,-0.5,0.0,0.5,0.9};
  double statevals[5]={0.05, 0.25, 0.5, 0.75, 0.95};

  int s0c,s1c,s2c,s3c;

  int score;

  thecart->nmarkov_long=false;
  thecart->generalization_test=false;

  //Evaluate each organism on a test
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {

    //shouldn't happen
    if (((*curorg)->gnome)==0) {
      cout<<"ERROR EMPTY GEMOME!"<<endl;
      cin>>pause;
    }

    if (pole2_evaluate((*curorg),velocity,thecart)) win=true;

  }

  //Average and max their fitnesses for dumping to file and snapshot
  for(curspecies=(pop->species).begin();curspecies!=(pop->species).end();++curspecies) {

    //This experiment control routine issues commands to collect ave
    //and max fitness, as opposed to having the snapshot do it, 
    //because this allows flexibility in terms of what time
    //to observe fitnesses at

    (*curspecies)->compute_average_fitness();
    (*curspecies)->compute_max_fitness();
  }

  //Take a snapshot of the population, so that it can be
  //visualized later on
  //if ((generation%1)==0)
  //  pop->snapshot();

  //Find the champion in the markov case simply for stat collection purposes
  if (thecart->MARKOV) {
    champ_fitness=0.0;
    for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg) {
      if (((*curorg)->fitness1)>champ_fitness) {
	champ=(*curorg);
	champ_fitness=champ->fitness1;
	champgenes=champ->gnome->genes.size();
	champnodes=champ->gnome->nodes.size();
	winnernum=champ->gnome->genome_id;
      }
    }
  }

  //Check for winner in Non-Markov case
  if (!(thecart->MARKOV)) {
    
    cout<<"Non-markov case"<<endl;

    //Sort the species
    for(curspecies=(pop->species).begin();curspecies!=(pop->species).end();++curspecies) {
      sorted_species.push_back(*curspecies);
    }

    //sorted_species.sort(order_new_species);
    std::sort(sorted_species.begin(), sorted_species.end(), NEAT::order_new_species);

    //std::sort(sorted_species.begin(), sorted_species.end(), order_species);


    cout<<"Number of species sorted: "<<sorted_species.size()<<endl;

    //First update what is checked and unchecked
    for(curspecies=sorted_species.begin();curspecies!=sorted_species.end();++curspecies) {
      if (((*curspecies)->compute_max_fitness())>((*curspecies)->max_fitness_ever))
	(*curspecies)->checked=false;

    }

    //Now find a species that is unchecked
    curspecies=sorted_species.begin();
    cout<<"Is the first species checked? "<<(*curspecies)->checked<<endl;
    while((curspecies!=(sorted_species.end()))&&
	  ((*curspecies)->checked))
    {
      cout<<"Species #"<<(*curspecies)->id<<" is checked"<<endl;
      ++curspecies;
    }

    if (curspecies==(sorted_species.end())) curspecies=sorted_species.begin();

    //Remember it was checked
    (*curspecies)->checked=true;
    cout<<"Is the species now checked? "<<(*curspecies)->checked<<endl;

    //Extract the champ
    cout<<"Champ chosen from Species "<<(*curspecies)->id<<endl;
    champ=(*curspecies)->get_champ();
    champ_fitness=champ->fitness1;
    cout<<"Champ is organism #"<<champ->gnome->genome_id<<endl;
    cout<<"Champ fitness: "<<champ_fitness<<endl;
    winnernum=champ->gnome->genome_id;

    cout<<champ->gnome<<endl;

    //Now check to make sure the champ can do 100,000
    thecart->nmarkov_long=true;
    thecart->generalization_test=false;

    //The champ needs tp be flushed here because it may have
    //leftover activation from its last test run that could affect
    //its recurrent memory
    (champ->net)->flush();


    //champ->gnome->print_to_filename("tested");
    
    if (pole2_evaluate(champ,velocity,thecart)) {
      cout<<"The champ passed the 100,000 test!"<<endl;

      thecart->nmarkov_long=false;

      //Given that the champ passed, now run it on generalization tests
      score=0;
      for (s0c=0;s0c<=4;++s0c)
	for (s1c=0;s1c<=4;++s1c)
	  for (s2c=0;s2c<=4;++s2c)
	    for (s3c=0;s3c<=4;++s3c) {
	      thecart->state[0] = statevals[s0c] * 4.32 - 2.16;
	      thecart->state[1] = statevals[s1c] * 2.70 - 1.35;
	      thecart->state[2] = statevals[s2c] * 0.12566304 - 0.06283152;
	      /* 0.06283152 =  3.6 degrees */
	      thecart->state[3] = statevals[s3c] * 0.30019504 - 0.15009752;
	      /* 00.15009752 =  8.6 degrees */
	      thecart->state[4]=0.0;
	      thecart->state[5]=0.0;
	      
	      cout<<"On combo "<<thecart->state[0]<<" "<<thecart->state[1]<<" "<<thecart->state[2]<<" "<<thecart->state[3]<<endl;
	      thecart->generalization_test=true;
	      
	      (champ->net)->flush();  //Reset the champ for each eval

	      if (pole2_evaluate(champ,velocity,thecart)) {
		cout<<"----------------------------The champ passed its "<<score<<"th test"<<endl;
		score++;
	      }
	      
	    }

      if (score>=200) {
	cout<<"The champ wins!!! (generalization = "<<score<<" )"<<endl;
	oFile<<"(generalization = "<<score<<" )"<<endl;
	oFile<<"generation= "<<generation<<endl;
        (champ->gnome)->print_to_file(oFile);
	champ_fitness=champ->fitness1;
	champgenes=champ->gnome->genes.size();
	champnodes=champ->gnome->nodes.size();
	winnernum=champ->gnome->genome_id;
	win=true;
      }
      else {
	cout<<"The champ couldn't generalize"<<endl;
	champ->fitness1=champ_fitness; //Restore the champ's fitness
      }
    }
    else {
      cout<<"The champ failed the 100,000 test :("<<endl;
      cout<<"made score "<<champ->fitness1<<endl;
      champ->fitness1=champ_fitness; //Restore the champ's fitness
    }
  }
  
  //Only print to file every print_every generations
  if  (win||
       ((generation%(NEAT::print_every))==0)) {
    cout<<"printing file: "<<filename<<endl;
    pop->print_to_file_by_species(filename);
  }

  if ((win)&&((pop->winnergen)==0)) pop->winnergen=generation;

  //Prints a champion out on each generation
  //IMPORTANT: This causes generational file output!
  print_Genome_tofile(champ->gnome,"champ");

  //Create the next generation
  pop->epoch(generation, filename);

  return (int) champ_fitness;
}

bool pole2_evaluate(Organism *org,bool velocity, CartPole *thecart) {
  Network *net;

  int thresh;  /* How many visits will be allowed before giving up 
		  (for loop detection)  NOW OBSOLETE */

  int pause;

  net=org->net;

  thresh=100;  //this is obsolete

  //DEBUG :  Check flushedness of org
  //org->net->flush_check();

  //Try to balance a pole now
  org->fitness1 = thecart->evalNet(net,thresh);

#ifndef NO_SCREEN_OUT
  if (org->pop_champ_child)
    cout<<" <<DUPLICATE OF CHAMPION>> ";

  //Output to screen
  cout<<"Org "<<(org->gnome)->genome_id<<" fitness: "<<org->fitness1;
  cout<<" ("<<(org->gnome)->genes.size();
  cout<<" / "<<(org->gnome)->nodes.size()<<")";
  cout<<"   ";
  if (org->mut_struct_baby) cout<<" [struct]";
  if (org->mate_baby) cout<<" [mate]";
  cout<<endl;
#endif

  if ((!(thecart->generalization_test))&&(!(thecart->nmarkov_long)))
  if (org->pop_champ_child) {
    cout<<org->gnome<<endl;
    //DEBUG CHECK
    if (org->high_fit>org->fitness1) {
      cout<<"ALERT: ORGANISM DAMAGED"<<endl;
      print_Genome_tofile(org->gnome,"failure_champ_genome");
      cin>>pause;
    }
  }

  //Decide if its a winner, in Markov Case
  if (thecart->MARKOV) {
    if (org->fitness1>=(thecart->maxFitness-1)) { 
      org->winner=true;
      return true;
    }
    else {
      org->winner=false;
      return false;
    }
  }
  //if doing the long test non-markov 
  else if (thecart->nmarkov_long) {
    if (org->fitness1>=99999) { 
      //if (org->fitness1>=9000) { 
      org->winner=true;
      return true;
    }
    else {
      org->winner=false;
      return false;
    }
  }
  else if (thecart->generalization_test) {
    if (org->fitness1>=999) {
      org->winner=true;
      return true;
    }
    else {
      org->winner=false;
      return false;
    }
  }
  else {
    org->winner=false;
    return false;  //Winners not decided here in non-Markov
  }
}

CartPole::CartPole(bool randomize,bool velocity)
{
  maxFitness = 100000;

  MARKOV=velocity;

  MIN_INC = 0.001;
  POLE_INC = 0.05;
  MASS_INC = 0.01;

  LENGTH_2 = 0.05;
  MASSPOLE_2 = 0.01;

  // CartPole::reset() which is called here
}

//Faustino Gomez wrote this physics code using the differential equations from 
//Alexis Weiland's paper and added the Runge-Kutta himself.
double CartPole::evalNet(Network *net,int thresh)
{
  int steps=0;
  double input[NUM_INPUTS];
  double output;

  int nmarkovmax;  

  double nmarkov_fitness;

  double jiggletotal; //total jiggle in last 100
  int count;  //step counter

  //init(randomize);		// restart at some point
  
  if (nmarkov_long) nmarkovmax=100000;
  else if (generalization_test) nmarkovmax=1000;
  else nmarkovmax=1000;


  init(0);

  if (MARKOV) {
    while (steps++ < maxFitness) {
      
         
      input[0] = state[0] / 4.8;
      input[1] = state[1] /2;
      input[2] = state[2]  / 0.52;
      input[3] = state[3] /2;
      input[4] = state[4] / 0.52;
      input[5] = state[5] /2;
      input[6] = .5;
      
      net->load_sensors(input);
      
      //Activate the net
      //If it loops, exit returning only fitness of 1 step
      if (!(net->activate())) return 1.0;
      
      output=(*(net->outputs.begin()))->activation;
      
      performAction(output,steps);
      
      if (outsideBounds())	// if failure
	break;			// stop it now
    }
    return (double) steps;
  }
  else {  //NON MARKOV CASE

    while (steps++ < nmarkovmax) {
      

     //Do special parameter summing on last hundred
     //if ((steps==900)&&(!nmarkov_long)) last_hundred=true;

     /*
     input[0] = state[0] / 4.8;
     input[1] = 0.0;
     input[2] = state[2]  / 0.52;
     input[3] = 0.0;
     input[4] = state[4] / 0.52;
     input[5] = 0.0;
     input[6] = .5;
     */

      //cout<<"nmarkov_long: "<<nmarkov_long<<endl;

      //if (nmarkov_long)
      //cout<<"step: "<<steps<<endl;

     input[0] = state[0] / 4.8;
     input[1] = state[2]  / 0.52;
     input[2] = state[4] / 0.52;
     input[3] = .5;
      
      net->load_sensors(input);

      //cout<<"inputs: "<<input[0]<<" "<<input[1]<<" "<<input[2]<<" "<<input[3]<<endl;

      //Activate the net
      //If it loops, exit returning only fitness of 1 step
      if (!(net->activate())) return 0.0001;
      
      output=(*(net->outputs.begin()))->activation;

      //cout<<"output: "<<output<<endl;

      performAction(output,steps);

      if (outsideBounds())	// if failure
	break;			// stop it now

      if (nmarkov_long&&(outsideBounds()))	// if failure
	break;			// stop it now
    }

   //If we are generalizing we just need to balance it a while
   if (generalization_test)
     return (double) balanced_sum;
 
   //Sum last 100
   if ((steps>100)&&(!nmarkov_long)) {

     jiggletotal=0;
     cout<<"step "<<steps-99-2<<" to step "<<steps-2<<endl;
     //Adjust for array bounds and count
     for (count=steps-99-2;count<=steps-2;count++)
       jiggletotal+=jigglestep[count];
   }

   if (!nmarkov_long) {
     if (balanced_sum>100) 
       nmarkov_fitness=((0.1*(((double) balanced_sum)/1000.0))+
			(0.9*(0.75/(jiggletotal))));
     else nmarkov_fitness=(0.1*(((double) balanced_sum)/1000.0));

#ifndef NO_SCREEN_OUTR
     cout<<"Balanced:  "<<balanced_sum<<" jiggle: "<<jiggletotal<<" ***"<<endl;
#endif

     return nmarkov_fitness;
   }
   else return (double) steps;

  }

}

void CartPole::init(bool randomize)
{
  static int first_time = 1;

  if (!MARKOV) {
    //Clear all fitness records
    cartpos_sum=0.0;
    cartv_sum=0.0;
    polepos_sum=0.0;
    polev_sum=0.0;
  }

  balanced_sum=0; //Always count # balanced

  last_hundred=false;

  /*if (randomize) {
    state[0] = (lrand48()%4800)/1000.0 - 2.4;
    state[1] = (lrand48()%2000)/1000.0 - 1;
    state[2] = (lrand48()%400)/1000.0 - 0.2;
    state[3] = (lrand48()%400)/1000.0 - 0.2;
    state[4] = (lrand48()%3000)/1000.0 - 1.5;
    state[5] = (lrand48()%3000)/1000.0 - 1.5;
  }
  else {*/


  if (!generalization_test) {
    state[0] = state[1] = state[3] = state[4] = state[5] = 0;
    state[2] = 0.07; // one_degree;
  }
  else {
    state[4] = state[5] = 0;
  }

    //}
  if(first_time){
    cout<<"Initial Long pole angle = %f\n"<<state[2]<<endl;;
    cout<<"Initial Short pole length = %f\n"<<LENGTH_2<<endl;
    first_time = 0;
  }
}

void CartPole::performAction(double output, int stepnum)
{ 
  
  int i;
  double  dydx[6];

  const bool RK4=true; //Set to Runge-Kutta 4th order integration method
  const double EULER_TAU= TAU/4;
 
  /*random start state for long pole*/
  /*state[2]= drand48();   */
     
  /*--- Apply action to the simulated cart-pole ---*/

  if(RK4){
    for(i=0;i<2;++i){
      dydx[0] = state[1];
      dydx[2] = state[3];
      dydx[4] = state[5];
      step(output,state,dydx);
      rk4(output,state,dydx,state);
    }
  }
  else{
    for(i=0;i<8;++i){
      step(output,state,dydx);
      state[0] += EULER_TAU * dydx[0];
      state[1] += EULER_TAU * dydx[1];
      state[2] += EULER_TAU * dydx[2];
      state[3] += EULER_TAU * dydx[3];
      state[4] += EULER_TAU * dydx[4];
      state[5] += EULER_TAU * dydx[5];
    }
  }

  //Record this state
  cartpos_sum+=fabs(state[0]);
  cartv_sum+=fabs(state[1]);
  polepos_sum+=fabs(state[2]);
  polev_sum+=fabs(state[3]);
  if (stepnum<=1000)
    jigglestep[stepnum-1]=fabs(state[0])+fabs(state[1])+fabs(state[2])+fabs(state[3]);

  if (false) {
    //cout<<"[ x: "<<state[0]<<" xv: "<<state[1]<<" t1: "<<state[2]<<" t1v: "<<state[3]<<" t2:"<<state[4]<<" t2v: "<<state[5]<<" ] "<<
    //cartpos_sum+cartv_sum+polepos_sum+polepos_sum+polev_sum<<endl;
    if (!(outsideBounds())) {
      if (balanced_sum<1000) {
	cout<<".";
	++balanced_sum;
      }
    }
    else {
      if (balanced_sum==1000)
	balanced_sum=1000;
      else balanced_sum=0;
    }
  }
  else if (!(outsideBounds()))
    ++balanced_sum;

}

void CartPole::step(double action, double *st, double *derivs)
{
    double force,costheta_1,costheta_2,sintheta_1,sintheta_2,
          gsintheta_1,gsintheta_2,temp_1,temp_2,ml_1,ml_2,fi_1,fi_2,mi_1,mi_2;

    force =  (action - 0.5) * FORCE_MAG * 2;
    costheta_1 = cos(st[2]);
    sintheta_1 = sin(st[2]);
    gsintheta_1 = GRAVITY * sintheta_1;
    costheta_2 = cos(st[4]);
    sintheta_2 = sin(st[4]);
    gsintheta_2 = GRAVITY * sintheta_2;
    
    ml_1 = LENGTH_1 * MASSPOLE_1;
    ml_2 = LENGTH_2 * MASSPOLE_2;
    temp_1 = MUP * st[3] / ml_1;
    temp_2 = MUP * st[5] / ml_2;
    fi_1 = (ml_1 * st[3] * st[3] * sintheta_1) +
           (0.75 * MASSPOLE_1 * costheta_1 * (temp_1 + gsintheta_1));
    fi_2 = (ml_2 * st[5] * st[5] * sintheta_2) +
           (0.75 * MASSPOLE_2 * costheta_2 * (temp_2 + gsintheta_2));
    mi_1 = MASSPOLE_1 * (1 - (0.75 * costheta_1 * costheta_1));
    mi_2 = MASSPOLE_2 * (1 - (0.75 * costheta_2 * costheta_2));
    
    derivs[1] = (force + fi_1 + fi_2)
                 / (mi_1 + mi_2 + MASSCART);
    
    derivs[3] = -0.75 * (derivs[1] * costheta_1 + gsintheta_1 + temp_1)
                 / LENGTH_1;
    derivs[5] = -0.75 * (derivs[1] * costheta_2 + gsintheta_2 + temp_2)
                  / LENGTH_2;

}

void CartPole::rk4(double f, double y[], double dydx[], double yout[])
{

	int i;

	double hh,h6,dym[6],dyt[6],yt[6];


	hh=TAU*0.5;
	h6=TAU/6.0;
	for (i=0;i<=5;i++) yt[i]=y[i]+hh*dydx[i];
	step(f,yt,dyt);
	dyt[0] = yt[1];
	dyt[2] = yt[3];
	dyt[4] = yt[5];
	for (i=0;i<=5;i++) yt[i]=y[i]+hh*dyt[i];
	step(f,yt,dym);
	dym[0] = yt[1];
	dym[2] = yt[3];
	dym[4] = yt[5];
	for (i=0;i<=5;i++) {
		yt[i]=y[i]+TAU*dym[i];
		dym[i] += dyt[i];
	}
	step(f,yt,dyt);
	dyt[0] = yt[1];
	dyt[2] = yt[3];
	dyt[4] = yt[5];
	for (i=0;i<=5;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
}

bool CartPole::outsideBounds()
{
  const double failureAngle = thirty_six_degrees; 

  return 
    state[0] < -2.4              || 
    state[0] > 2.4               || 
    state[2] < -failureAngle     ||
    state[2] > failureAngle      ||
    state[4] < -failureAngle     ||
    state[4] > failureAngle;  
}

void CartPole::nextTask()
{

   LENGTH_2 += POLE_INC;   /* LENGTH_2 * INCREASE;   */
   MASSPOLE_2 += MASS_INC; /* MASSPOLE_2 * INCREASE; */
   //  ++new_task;
   cout<<"#Pole Length %2.4f\n"<<LENGTH_2<<endl;
}

void CartPole::simplifyTask()
{
  if(POLE_INC > MIN_INC) {
    POLE_INC = POLE_INC/2;
    MASS_INC = MASS_INC/2;
    LENGTH_2 -= POLE_INC;
    MASSPOLE_2 -= MASS_INC;
    cout<<"#SIMPLIFY\n"<<endl;
    cout<<"#Pole Length %2.4f\n"<<LENGTH_2;
  }
  else
    {
      cout<<"#NO TASK CHANGE\n"<<endl;
    }
}
