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
#include "network.h"
#include <iostream>
#include <sstream>

using namespace NEAT;

Network::Network(std::vector<NNode*> in,std::vector<NNode*> out,std::vector<NNode*> all,int netid) {
  inputs=in;
  outputs=out;
  all_nodes=all;
  name=0;   //Defaults to no name  ..NOTE: TRYING TO PRINT AN EMPTY NAME CAN CAUSE A CRASH
  numnodes=-1;
  numlinks=-1;
  net_id=netid;
  adaptable=false;
}

Network::Network(std::vector<NNode*> in,std::vector<NNode*> out,std::vector<NNode*> all,int netid, bool adaptval) {
  inputs=in;
  outputs=out;
  all_nodes=all;
  name=0;   //Defaults to no name  ..NOTE: TRYING TO PRINT AN EMPTY NAME CAN CAUSE A CRASH                                    
  numnodes=-1;
  numlinks=-1;
  net_id=netid;
  adaptable=adaptval;
}


Network::Network(int netid) {
			name=0; //Defaults to no name
			numnodes=-1;
			numlinks=-1;
			net_id=netid;
			adaptable=false;
		}

Network::Network(int netid, bool adaptval) {
  name=0; //Defaults to no name                                                                                               
  numnodes=-1;
  numlinks=-1;
  net_id=netid;
  adaptable=adaptval;
}


Network::Network(const Network& network)
{
	std::vector<NNode*>::const_iterator curnode;

	// Copy all the inputs
	for(curnode = network.inputs.begin(); curnode != network.inputs.end(); ++curnode) {
		NNode* n = new NNode(**curnode);
		inputs.push_back(n);
		all_nodes.push_back(n);
	}

	// Copy all the outputs
	for(curnode = network.outputs.begin(); curnode != network.outputs.end(); ++curnode) {
		NNode* n = new NNode(**curnode);
		outputs.push_back(n);
		all_nodes.push_back(n);
	}

	if(network.name)
		name = strdup(network.name);
	else
		name = 0;

	numnodes = network.numnodes;
	numlinks = network.numlinks;
	net_id = network.net_id;
	adaptable = network.adaptable;
}

Network::~Network() {
			if (name!=0)
				delete [] name;

			destroy();  // Kill off all the nodes and links

		}

// Puts the network back into an initial state
void Network::flush() {
	std::vector<NNode*>::iterator curnode;

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
		(*curnode)->flushback();
	}
}

// Debugger: Checks network state
void Network::flush_check() {
	std::vector<NNode*>::iterator curnode;
	std::vector<NNode*>::iterator location;
	std::vector<NNode*> seenlist;  //List of nodes not to doublecount

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {    
        location= std::find(seenlist.begin(),seenlist.end(),(*curnode));
		if (location==seenlist.end()) {
			seenlist.push_back(*curnode);
			(*curnode)->flushback_check(seenlist);
		}
	}
}

// Aditya: Turn off all the outputs at the start of every activate (For RECALL task)
// Ensures that the output node is activated again with the incoming active values each time 
void Network::switch_outputsoff() {
	std::vector<NNode*>::iterator curnode;

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
		((*curnode)->activation_count)=0;
	}
}

// If all output are not active then return true
bool Network::outputsoff() {
	std::vector<NNode*>::iterator curnode;

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
		if (((*curnode)->activation_count)==0) return true;
	}

	return false;
}

// Print the connections weights to a file separated by only carriage returns
void Network::print_links_tofile(char *filename) {
	std::vector<NNode*>::iterator curnode;
	std::vector<Link*>::iterator curlink;

    std::ofstream oFile(filename);

	//Make sure it worked
	//if (!oFile) {
	//	cerr<<"Can't open "<<filename<<" for output"<<endl;
		//return 0;
	//}

	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
		if (((*curnode)->type)!=SENSOR) {
			for(curlink=((*curnode)->incoming).begin(); curlink!=((*curnode)->incoming).end(); ++curlink) {
                oFile << (*curlink)->in_node->node_id << " -> " <<( *curlink)->out_node->node_id << " : " << (*curlink)->weight << std::endl;
			} // end for loop on links
		} //end if
	} //end for loop on nodes

	oFile.close();

} //print_links_tofile

void Network::lstm_activate(NNode* curnode){
	//if (curnode->active_in_flag) {
                //std::cout<<"Activating "<<curnode->node_id<<" with "<<curnode->activesum<<": ";

		//Keep a memory of activations for potential time delayed connections
		curnode->last_activation2=curnode->last_activation;
		curnode->last_activation=curnode->activation;

                //Activate gate controls
                double rd_gate_activation =  NEAT::fsigmoid(curnode->activesum_rd,1.0,2.4621365);//Gate control values range 0-1(Can do binary)
                double wr_gate_activation =  NEAT::fsigmoid(curnode->activesum_wr,1.0,2.4621365);//Gate control values range 0-1(Can do binary)
                double fg_gate_activation =  NEAT::fsigmoid(curnode->activesum_fg,1.0,2.4621365);//Gate control values range 0-1(Can do binary)
                
                //double fg_gate_activation = 0.55;
                //Making control signals binary 
                if (rd_gate_activation>0.5) rd_gate_activation=1.0;
                else rd_gate_activation=0.0;
                if (wr_gate_activation>0.5) wr_gate_activation=1.0;
                else wr_gate_activation=0.0;
                //if (fg_gate_activation>0.5) fg_gate_activation=1.0;
                //else fg_gate_activation=0.0;
                //std::cout<<"rd: "<<curnode->activesum_rd<< " "<<rd_gate_activation<<std::endl;
                //std::cout<<"wr: "<<curnode->activesum_wr<< " "<<wr_gate_activation<<std::endl;
                //std::cout<<"fg: "<<curnode->activesum_fg<< " "<<fg_gate_activation<<std::endl;
                //std::cout<<"last_Cell_state: "<< curnode->lstm_cell_state<<std::endl;
                //std::cout<<"activesum: "<< curnode->activesum<<std::endl;
                //curnode->lstm_cell_state = (curnode->activesum)*wr_gate_activation + (curnode->lstm_cell_state)*fg_gate_activation; //Input and history gating
                if (curnode->lstm_cell_state == 0) {//ALlow new writes only when the LSTM cell is empty (write occurs only once for now)
                        curnode->lstm_cell_state = (curnode->activesum)*wr_gate_activation;// + (curnode->lstm_cell_state)*fg_gate_activation; //Input and history gating

                }
                else {//Gradually decay the value
                        curnode->lstm_cell_state = (curnode->lstm_cell_state)*fg_gate_activation; //Input and history gating
                } 
                //if (wr_gate_activation==1.0) {
                //        curnode->lstm_cell_state = (curnode->activesum);//If incoming write, overwrite the stored value
                //}
                //else {
                //        curnode->lstm_cell_state = curnode->lstm_cell_state;//If no incoming write, then maintain previously stored value
                //}
                //std::cout<<"new_Cell_state: "<< curnode->lstm_cell_state<<std::endl;

                double lstm_out =  NEAT::ftanh((curnode->lstm_cell_state),1.0,2.4621365);//(curnode->lstm_cell_state);//Always Read (No output gating) * rd_gate_activation; //Output gating
                //std::cout<<"lstm_out: "<< lstm_out<<std::endl;

                curnode->activation=lstm_out;//NEAT::ftanh(lstm_out,1.0,2.4621365); //Evolino uses tanh. Can try sigmoid as well
                //std::cout<<"lstm_final_out: "<< curnode->activation<<std::endl;

		//Increment the activation_count
		//First activation cannot be from nothing!!
		curnode->activation_count++;
                //curnode->active_out_flag = true;
	//}

}
void Network::setup_lstm_activate(NNode* curnode){
	curnode->activesum=0;
	curnode->activesum_rd=0;
	curnode->activesum_wr=0;
	curnode->activesum_fg=0;
	//curnode->active_in_flag=false;  //This will tell us if it has any active inputs
	std::vector<Link*>::iterator curlink;
        double add_amount = 0.0;
        //bool active_flag_rd=false, active_flag_wr=false, active_flag_fg=false, active_flag_inputdata=false;
	// For each incoming connection, add the activity from the connection to the activesum 
	for(curlink=(curnode->incoming).begin();curlink!=(curnode->incoming).end();++curlink) {
                add_amount = 0.0;
                add_amount=((*curlink)->weight)*(((*curlink)->in_node)->get_active_out());
		//std::cout<<"LSTM Node "<<(curnode)->node_id<<" adding "<<add_amount<<" from node "<<((*curlink)->in_node)->node_id<<std::endl;
                //std::cout<<"Link Type: "<<(*curlink)->link_gtype <<std::endl;
                if((*curlink)->link_gtype==NONE) {
			//if ((((*curlink)->in_node)->active_out_flag)) {
			        curnode->activesum+=add_amount;
                                //active_flag_inputdata=true;
		                //std::cout<<"LSTM Node Input Data activesum: "<<curnode->activesum<<" active_flag: "<<active_flag_inputdata <<std::endl;
                        //}
                }
                else if ((*curlink)->link_gtype==READ) {
			//if ((((*curlink)->in_node)->active_out_flag)) {
			        curnode->activesum_rd+=add_amount;
                                //active_flag_rd=true;
		                //std::cout<<"LSTM Node READ activesum: "<<curnode->activesum_rd<<" active_flag: "<<active_flag_rd <<std::endl;
                        //}
                }
                else if ((*curlink)->link_gtype==WRITE) {
			//if ((((*curlink)->in_node)->active_out_flag)) {
			        curnode->activesum_wr+=add_amount;
                                //active_flag_wr=true;
                                //std::cout<<"LSTM Node WRITE activesum: "<<curnode->activesum_wr<<" active_flag: "<<active_flag_wr <<std::endl;
                        //}
                }
                else if ((*curlink)->link_gtype==FORGET) {
			//if ((((*curlink)->in_node)->active_out_flag)) {
			        curnode->activesum_fg+=add_amount;
                                //active_flag_fg=true;
		                //std::cout<<"LSTM Node FORGET activesum: "<<curnode->activesum_fg<<" active_flag: "<<active_flag_fg <<std::endl;
                        //}
                }
        }
        //if (active_flag_inputdata) {//LSTM node is active only when write control is enabled 
                //curnode->active_in_flag=true;
        //}
                //std::cout<<" LSTM active_in_flag: "<<curnode->active_in_flag <<std::endl;
}
// Activates the net such that all outputs are active
// Returns true on success;
bool Network::activate() {
	std::vector<NNode*>::iterator curnode;
	std::vector<Link*>::iterator curlink;
	double add_amount;  //For adding to the activesum
	bool onetime; //Make sure we at least activate once
	int abortcount=0;  //Used in case the output is somehow truncated from the network

	//cout<<"Activating network: "<<this->genotype<<endl;

	//Keep activating until all the outputs have become active 
	//(This only happens on the first activation, because after that they
	// are always active)

        switch_outputsoff(); // Aditya: Ensures that the output node is activated again with the incoming active values each time 
	onetime=false;

	while(!onetime) {

		//++abortcount;

		//if (abortcount==20) {
		//	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
		//		//std::cout <<(*curnode)->node_id<< " :: "<<(*curnode)->active_flag ;  //This will tell us if it has any active inputs
		//		if ((int)(*curnode)->incoming.size() == 0) {
		//			//std::cout << " NO incoming connections";

		//		}
		//		//std::cout << std::endl;
		//	}
		//	//std::cout << " Net not activating"<<std::endl;
		//	return false;
		//	//cout<<"Inputs disconnected from output!"<<endl;
		//}
		//std::cout<<"Outputs are off"<<std::endl;

		// For each node, compute the sum of its incoming activation 
		for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
			//Ignore SENSORS

                        //std::cout<<"On node "<<(*curnode)->node_id<<std::endl;

			if (((*curnode)->type)!=SENSOR) {
                                if (((*curnode)->type)!=LSTM) {//For regular hidden neurons and output nodes
				        (*curnode)->activesum=0;
				        //(*curnode)->active_in_flag=false;  //This will tell us if it has any active inputs

				        // For each incoming connection, add the activity from the connection to the activesum 
				        for(curlink=((*curnode)->incoming).begin();curlink!=((*curnode)->incoming).end();++curlink) {
				        	//Handle possible time delays
				        	if (!((*curlink)->time_delay)) {
				        		add_amount=((*curlink)->weight)*(((*curlink)->in_node)->get_active_out());
				        		//if ((((*curlink)->in_node)->active_out_flag)){//||(((*curlink)->in_node)->type==SENSOR)) 
                                                                //(*curnode)->active_in_flag=true;
				        		        (*curnode)->activesum+=add_amount;
				        		        //std::cout<<"Node "<<(*curnode)->node_id<<" adding "<<add_amount<<" from node "<<((*curlink)->in_node)->node_id<<" active_in_flag: "<<(*curnode)->active_in_flag <<" "<<(*curnode)->activesum<<std::endl;
                                                        //}
				        	}
				        	else {
				        		//Input over a time delayed connection
				        		add_amount=((*curlink)->weight)*(((*curlink)->in_node)->get_active_out_td());
				        		(*curnode)->activesum+=add_amount;
				        	}

				        } //End for over incoming links
                                }//End if !LSTM
                                else {//For LSTM nodes
                                        setup_lstm_activate((*curnode));

                                }
			} //End if (((*curnode)->type)!=SENSOR) 

		} //End for over all nodes

		// Now activate all the non-sensor nodes off their incoming activation 
		for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {

                                //std::cout<<"Node type: "<<((*curnode)->type)<<std::endl;
			if (((*curnode)->type)!=SENSOR) {
                                if (((*curnode)->type)!=LSTM) {//For regular hidden neurons and output nodes
				        //Only activate if some active input came in
				        //if ((*curnode)->active_in_flag) {
                                                //std::cout<<"Activating "<<(*curnode)->node_id<<" with "<<(*curnode)->activesum<<": ";

				        	//Keep a memory of activations for potential time delayed connections
				        	(*curnode)->last_activation2=(*curnode)->last_activation;
				        	(*curnode)->last_activation=(*curnode)->activation;

				        	//If the node is being overrided from outside,
				        	//stick in the override value
				        	if ((*curnode)->overridden()) {
				        		//Set activation to the override value and turn off override
				        		(*curnode)->activate_override();
				        	}
				        	else {
				        		//Now run the net activation through an activation function
				        		if ((*curnode)->ftype==SIGMOID) {
                                                                (*curnode)->activation=NEAT::fReLu((*curnode)->activesum,1.0,2.4621365);  //Sigmoidal activation- see comments under fsigmoid //Changed slope from 4.924273 to 1.0 to allow for 4 different output node categories
                                                        }
                                                        else {
                                                                (*curnode)->activation=NEAT::fReLu((*curnode)->activesum,1.0,2.4621365);  //Rectified Linear Units {max(0,x)}
                                                        }
                                                }
                                                //std::cout<<(*curnode)->activation<<std::endl;

				        	//Increment the activation_count
				        	//First activation cannot be from nothing!!
				        	(*curnode)->activation_count++;
                                                //(*curnode)->active_out_flag = true;
				        //}
                                        //else {
                                                //(*curnode)->active_out_flag = false;

                                        //}
                                }
                                else { //For LSTM Nodes
                                        lstm_activate((*curnode));
                                }
			}
                        //Aditya: For the first time activate is called after network flush, this
                        //prevents output getting activated twice with the same input value
                        //Sensor active flag is set every time sensors are loaded with new values
                        //else if (((*curnode)->type) == SENSOR && ((*curnode)->active_out_flag == true)) {
                                //std::cout<<" Disabling input node: "<<(*curnode)->node_id<<std::endl;
                                //(*curnode)->active_out_flag = false;
                        //}
		}

		onetime=true;
	}

	if (adaptable) {

	  //std::cout << "ADAPTING" << std:endl;

	  // ADAPTATION:  Adapt weights based on activations 
	  for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
	    //Ignore SENSORS
	    
	    //cout<<"On node "<<(*curnode)->node_id<<endl;
	    
	    if (((*curnode)->type)!=SENSOR) {
	      
	      // For each incoming connection, perform adaptation based on the trait of the connection 
	      for(curlink=((*curnode)->incoming).begin();curlink!=((*curnode)->incoming).end();++curlink) {
		
		if (((*curlink)->trait_id==2)||
		    ((*curlink)->trait_id==3)||
		    ((*curlink)->trait_id==4)) {
		  
		  //In the recurrent case we must take the last activation of the input for calculating hebbian changes
		  if ((*curlink)->is_recurrent) {
		    (*curlink)->weight=
		      hebbian((*curlink)->weight,maxweight,
			      (*curlink)->in_node->last_activation, 
			      (*curlink)->out_node->get_active_out(),
			      (*curlink)->params[0],(*curlink)->params[1],
			      (*curlink)->params[2]);
		    
		    
		  }
		  else { //non-recurrent case
		    (*curlink)->weight=
		      hebbian((*curlink)->weight,maxweight,
			      (*curlink)->in_node->get_active_out(), 
			      (*curlink)->out_node->get_active_out(),
			      (*curlink)->params[0],(*curlink)->params[1],
			      (*curlink)->params[2]);
		  }
		}
		
	      }
	      
	    }
	    
	  }
	  
	} //end if (adaptable)

	return true;  
}

//Recursive finds activation value of a node. Activates any node only once.
//Works only when there is no recurrent/loopy connection 
void Network::recursive_activation(NNode* curnode, bool &network_error){
	std::vector<Link*>::iterator curlink;
        //curnode-> active_flag = false;
        curnode-> activesum = 0.0;
        if (curnode->visited==false) {//Entering the node for the first time in this recursion (starting from a given output node)
                
                curnode->visited = true; //Entering the node recursion. This is set to prevent loop activations
                if ((curnode->type)==SENSOR) { //Sensor has already been activated by loading sensor values
                        curnode-> active_in_flag = true;
                }
                else if ((curnode-> active_in_flag)==false) {//If current node has not been activated before through a different path
	                curnode->activesum=0;
	                curnode->activesum_rd=0;
	                curnode->activesum_wr=0;
	                curnode->activesum_fg=0;
                        double add_amount = 0.0;
                        bool active_flag_rd, active_flag_wr, active_flag_fg, active_flag_inputdata=false;
                        for(curlink=(curnode->incoming).begin();curlink!=(curnode->incoming).end();++curlink) {
                                recursive_activation((*curlink)->in_node, network_error); //Recursive function
                                add_amount=((*curlink)->weight)*(((*curlink)->in_node)->get_active_out());
                                if ((*curlink)->in_node->active_in_flag) {//Only if the node is not floating, use its value 
                                        //curnode->activesum+= (((*curlink)->in_node)->get_active_out())*((*curlink)->weight);
                                        if((*curlink)->link_gtype==NONE) {
	                                	curnode->activesum+=add_amount;
	                                	active_flag_inputdata=true;
	                                        //std::cout<<"LSTM/Regular Node Input Data activesum: "<<curnode->activesum<<" active_flag: "<<active_flag_inputdata <<std::endl;
                                        }
                                        else if ((*curlink)->link_gtype==READ) {
	                                	curnode->activesum_rd+=add_amount;
	                                	active_flag_rd=true;
	                                        //std::cout<<"LSTM Node READ activesum: "<<curnode->activesum_rd<<" active_flag: "<<active_flag_rd <<std::endl;
                                        }
                                        else if ((*curlink)->link_gtype==WRITE) {
	                                	curnode->activesum_wr+=add_amount;
	                                	active_flag_wr=true;
	                                        //std::cout<<"LSTM Node WRITE activesum: "<<curnode->activesum_wr<<" active_flag: "<<active_flag_rd <<std::endl;
                                        }
                                        else if ((*curlink)->link_gtype==FORGET) {
	                                	curnode->activesum_fg+=add_amount;
	                                	active_flag_fg=true;
	                                        //std::cout<<"LSTM Node FORGET activesum: "<<curnode->activesum_fg<<" active_flag: "<<active_flag_fg <<std::endl;
                                        }
                                }
                        }
                        if ((curnode->type)!=LSTM) {//For Regular nodes, just check the incoming data 
                                if (active_flag_inputdata){
                                        curnode-> active_in_flag = true; //This node is active if at least one incoming link is active
                                }
                        }
                        else {//For LSTM nodes, make sure both the data and control are active
                                if (active_flag_rd && active_flag_wr && active_flag_fg && active_flag_inputdata){
                                        curnode-> active_in_flag = true; 
                                }
                        }
                        if ((curnode-> active_in_flag)) {
                                if ((curnode->type)!=LSTM) {//Activate Regular nodes
	        	                //Now run the net activation through an activation function
	        	                if (curnode->ftype==SIGMOID) {
                                                curnode->activation=NEAT::fsigmoid(curnode->activesum,1.0,2.4621365);  //Sigmoidal activation- see comments under fsigmoid //Changed slope from 4.924273 to 1.0 to allow for 4 different output node categories
                                        }
                                        else {
                                                curnode->activation=NEAT::fReLu(curnode->activesum,1.0,2.4621365);  //Rectified Linear Units {max(0,x)}
                                        }
	                                curnode->activation_count++;
                                }
                                else {//Activate LSTM node
                                        lstm_activate(curnode);
                                }
                                //Check to see if this node has been activated more than once
                                if (curnode->activation_count > 1) {
                                        print_links_tofile("error_node_multi_activations.txt");
                                        std::cout<< " ERRORR: Node: " <<curnode->node_id<<" being activated more than once"<<std::endl;
                                        network_error = true;
                                        //exit(0);
                                }
                        }
                        else {//Floating hidden nodes are OK. Just ignore them
                                //std::cout<<" WARNINGGGGGGGGGGGGGGGGGGGGGGG: Node: "<< curnode->node_id<<" cannot be floating. Make sure memory_startgenes has valid paths to output. Otherwise, could be a BUG (NOT NECESSARILY)"<<std::endl;
                                //exit(0);
                        }
                }
                curnode->visited = false; //Exiting the node
        }
}

// Activates the net such that all paths to the output are activated just once. 
// Used for static problems like image classification
// Works only when there is no recurrent/loopy connection 
// Returns true on success;
bool Network::activate_static(bool &network_error) {
	std::vector<NNode*>::iterator curnode;
	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
                recursive_activation(*curnode, network_error);
	}
        if (outputsoff()) {//If any output is floating
                return false;
        }
        else {
                return true;
        }
}
// THIS WAS NOT USED IN THE FINAL VERSION, AND NOT FULLY IMPLEMENTED,   
// BUT IT SHOWS HOW SOMETHING LIKE THIS COULD BE INITIATED
// Note that checking networks for loops in general in not necessary
// and therefore I stopped writing this function
// Check Network for loops.  Return true if its ok, false if there is a loop.
//bool Network::integrity() {
//  std::vector<NNode*>::iterator curnode;
//  std::vector<std::vector<NNode*>*> paths;
//  int count;
//  std::vector<NNode*> *newpath;
//  std::vector<std::vector<NNode*>*>::iterator curpath;

//  for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
//    newpath=new std::vector<NNode*>();
//    paths.push_back(newpath);
//    if (!((*curnode)->integrity(newpath))) return false;
//  }

//Delete the paths now that we are done
//  curpath=paths.begin();
//  for(count=0;count<paths.size();count++) {
//    delete (*curpath);
//    curpath++;
//  }

//  return true;
//}

// Prints the values of its outputs
void Network::show_activation() {
	std::vector<NNode*>::iterator curnode;
	int count;

	//if (name!=0)
	//  cout<<"Network "<<name<<" with id "<<net_id<<" outputs: (";
	//else cout<<"Network id "<<net_id<<" outputs: (";

	count=1;
	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
		//cout<<"[Output #"<<count<<": "<<(*curnode)<<"] ";
		count++;
	}

	//cout<<")"<<endl;
}

void Network::show_input() {
	std::vector<NNode*>::iterator curnode;
	int count;

	//if (name!=0)
	//  cout<<"Network "<<name<<" with id "<<net_id<<" inputs: (";
	//else cout<<"Network id "<<net_id<<" outputs: (";

	count=1;
	for(curnode=inputs.begin();curnode!=inputs.end();++curnode) {
		//cout<<"[Input #"<<count<<": "<<(*curnode)<<"] ";
		count++;
	}

	//cout<<")"<<endl;
}

// Add an input
void Network::add_input(NNode *in_node) {
	inputs.push_back(in_node);
}

// Add an output
void Network::add_output(NNode *out_node) {
	outputs.push_back(out_node);
}

// Takes an array of sensor values and loads it into SENSOR inputs ONLY
void Network::load_sensors(double *sensvals) {
	//int counter=0;  //counter to move through array
	std::vector<NNode*>::iterator sensPtr;

	for(sensPtr=inputs.begin();sensPtr!=inputs.end();++sensPtr) {
		//only load values into SENSORS (not BIASes)
		if (((*sensPtr)->type)==SENSOR) {
			(*sensPtr)->sensor_load(*sensvals);
			sensvals++;
		}
	}
}

void Network::load_sensors(const std::vector<double> &sensvals) {
	//int counter=0;  //counter to move through array
	std::vector<NNode*>::iterator sensPtr;
	std::vector<double>::const_iterator valPtr;

	for(valPtr = sensvals.begin(), sensPtr = inputs.begin(); sensPtr != inputs.end() && valPtr != sensvals.end(); ++sensPtr, ++valPtr) {
		//only load values into SENSORS (not BIASes)
		if (((*sensPtr)->type)==SENSOR) {
			(*sensPtr)->sensor_load(*valPtr);
			//sensvals++;
		}
	}
}


// Takes and array of output activations and OVERRIDES 
// the outputs' actual activations with these values (for adaptation)
void Network::override_outputs(double* outvals) {

	std::vector<NNode*>::iterator outPtr;

	for(outPtr=outputs.begin();outPtr!=outputs.end();++outPtr) {
		(*outPtr)->override_output(*outvals);
		outvals++;
	}

}

void Network::give_name(char *newname) {
	char *temp;
	char *temp2;
	temp=new char[strlen(newname)+1];
	strcpy(temp,newname);
	if (name==0) name=temp;
	else {
		temp2=name;
		delete temp2;
		name=temp;
	}
}

//Make a list of active nodes and links
void Network::find_active_paths() {
	int counter=0;
	std::vector<NNode*>::iterator curnode;
	
        //Reset variables
        for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
                (*curnode)->visited = false;
                (*curnode)->on_active_path = false;
        }

        //Recurse from output to input
	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
	        find_active_paths_helper((*curnode));
	}
}

bool Network::find_active_paths_helper(NNode *curnode) {
	std::vector<Link*> innodes=curnode->incoming;
	std::vector<Link*>::iterator curlink;
	std::vector<NNode*>::iterator location;
        bool temp = false;

        //If the node is sensor or already has been found to be on active path, 
        //then propagate this information up to the outputs
        if (((curnode->type)==SENSOR) || (curnode->visited == true && curnode->on_active_path==true)) { 
                curnode->on_active_path = true; //For sensor
                return true;
        }

        //Recurse deeper towards the input
        //If any branch is active, then this node/link are active as well
        else if (curnode->visited == false) {
                for(curlink=(curnode->incoming).begin();curlink!=(curnode->incoming).end();++curlink) {
                        temp = find_active_paths_helper((*curlink)->in_node);
                        if (temp == true) {//If any branch is active, then this node/link are active as well
                               curnode->on_active_path = true;
                               (*curlink)->on_active_path = true;
                        }
                }
                curnode->visited = true;
                if (curnode->on_active_path == true){
                        return true;
                }
                else {
                        return false;
                }
        }
}

// The following two methods recurse through a network from outputs
// down in order to count the number of nodes and links in the network.
// This can be useful for debugging genotype->phenotype spawning 
// (to make sure their counts correspond)

int Network::nodecount() {
	int counter=0;
	std::vector<NNode*>::iterator curnode;
	std::vector<NNode*>::iterator location;
	std::vector<NNode*> seenlist;  //List of nodes not to doublecount

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {

        location = std::find(seenlist.begin(),seenlist.end(),(*curnode));
		if (location==seenlist.end()) {
			counter++;
			seenlist.push_back(*curnode);
			nodecounthelper((*curnode),counter,seenlist);
		}
	}

	numnodes=counter;

	return counter;

}

void Network::nodecounthelper(NNode *curnode,int &counter,std::vector<NNode*> &seenlist) {
	std::vector<Link*> innodes=curnode->incoming;
	std::vector<Link*>::iterator curlink;
	std::vector<NNode*>::iterator location;

	if (!((curnode->type)==SENSOR)) {
		for(curlink=innodes.begin();curlink!=innodes.end();++curlink) {
            location= std::find(seenlist.begin(),seenlist.end(),((*curlink)->in_node));
			if (location==seenlist.end()) {
				counter++;
				seenlist.push_back((*curlink)->in_node);
				nodecounthelper((*curlink)->in_node,counter,seenlist);
			}
		}

	}

}

int Network::linkcount() {
	int counter=0;
	std::vector<NNode*>::iterator curnode;
	std::vector<NNode*> seenlist;  //List of nodes not to doublecount

	for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
		linkcounthelper((*curnode),counter,seenlist);
	}

	numlinks=counter;

	return counter;

}

void Network::linkcounthelper(NNode *curnode,int &counter,std::vector<NNode*> &seenlist) {
	std::vector<Link*> inlinks=curnode->incoming;
	std::vector<Link*>::iterator curlink;
	std::vector<NNode*>::iterator location;

    location = std::find(seenlist.begin(),seenlist.end(),curnode);
	if ((!((curnode->type)==SENSOR))&&(location==seenlist.end())) {
		seenlist.push_back(curnode);

		for(curlink=inlinks.begin();curlink!=inlinks.end();++curlink) {
			counter++;
			linkcounthelper((*curlink)->in_node,counter,seenlist);
		}

	}

}

// Destroy will find every node in the network and subsequently
// delete them one by one.  Since deleting a node deletes its incoming
// links, all nodes and links associated with a network will be destructed
// Note: Traits are parts of genomes and not networks, so they are not
//       deleted here
void Network::destroy() {
	std::vector<NNode*>::iterator curnode;
	std::vector<NNode*>::iterator location;
	std::vector<NNode*> seenlist;  //List of nodes not to doublecount

	// Erase all nodes from all_nodes list 

	for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode) {
		delete (*curnode);
	}


	// ----------------------------------- 

	//  OLD WAY-the old way collected the nodes together and then deleted them

	//for(curnode=outputs.begin();curnode!=outputs.end();++curnode) {
	//cout<<seenstd::vector<<endl;
	//cout<<curnode<<endl;
	//cout<<curnode->node_id<<endl;

	//  location=find(seenlist.begin(),seenlist.end(),(*curnode));
	//  if (location==seenlist.end()) {
	//    seenlist.push_back(*curnode);
	//    destroy_helper((*curnode),seenlist);
	//  }
	//}

	//Now destroy the seenlist, which is all the NNodes in the network
	//for(curnode=seenlist.begin();curnode!=seenlist.end();++curnode) {
	//  delete (*curnode);
	//}
}

void Network::destroy_helper(NNode *curnode,std::vector<NNode*> &seenlist) {
	std::vector<Link*> innodes=curnode->incoming;
	std::vector<Link*>::iterator curlink;
	std::vector<NNode*>::iterator location;

	if (!((curnode->type)==SENSOR)) {
		for(curlink=innodes.begin();curlink!=innodes.end();++curlink) {
            location = std::find(seenlist.begin(),seenlist.end(),((*curlink)->in_node));
			if (location==seenlist.end()) {
				seenlist.push_back((*curlink)->in_node);
				destroy_helper((*curlink)->in_node,seenlist);
			}
		}

	}

}

// This checks a POTENTIAL link between a potential in_node and potential out_node to see if it must be recurrent 
bool Network::is_recur(NNode *potin_node,NNode *potout_node,int &count,int thresh) {
	std::vector<Link*>::iterator curlink;


	++count;  //Count the node as visited

	if (count>thresh) {
		//cout<<"returning false"<<endl;
		return false;  //Short out the whole thing- loop detected
	}

	if (potin_node==potout_node) return true;
	else {
		//Check back on all links...
		for(curlink=(potin_node->incoming).begin();curlink!=(potin_node->incoming).end();curlink++) {
			//But skip links that are already recurrent
			//(We want to check back through the forward flow of signals only
			if (!((*curlink)->is_recurrent)) {
				if (is_recur((*curlink)->in_node,potout_node,count,thresh)) return true;
			}
		}
		return false;
	}
}

int Network::input_start() {
	input_iter=inputs.begin();
	return 1;
}

int Network::load_in(double d) {
	(*input_iter)->sensor_load(d);
	input_iter++;
	if (input_iter==inputs.end()) return 0;
	else return 1;
}


//Find the maximum number of neurons between an ouput and an input
int Network::max_depth() {
  std::vector<NNode*>::iterator curoutput; //The current output we are looking at
  int cur_depth; //The depth of the current node
  int max=0; //The max depth
  
  for(curoutput=outputs.begin();curoutput!=outputs.end();curoutput++) {
    cur_depth=(*curoutput)->depth(0,this);
    if (cur_depth>max) max=cur_depth;
  }

  return max;

}

