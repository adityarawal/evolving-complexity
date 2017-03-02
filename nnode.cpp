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
#include "nnode.h"
#include <iostream>
#include <sstream>
using namespace NEAT;

NNode::NNode(nodetype ntype,int nodeid) {
	active_in_flag=false;
	active_out_flag=false;
	visited=false;
	activesum=0;
	activesum_rd=0;
	activesum_wr=0;
	activesum_fg=0;
	activation=0;
        lstm_cell_state = 0.0;
	output=0;
	last_activation=0;
	last_activation2=0;
	type=ntype; //NEURON or SENSOR OR (LSTM NEURON) type
	activation_count=0; //Inactive upon creation
	node_id=nodeid;
	ftype=SIGMOID;
	nodetrait=0;
	gen_node_label=HIDDEN;
	dup=0;
	analogue=0;
	frozen=false;
        frozen_ip=false;//LSTM Input data. For LSTM, individual gates can be frozen separately
        frozen_rd=false;//LSTM Read gate. For LSTM, individual gates can be frozen separately
        frozen_wr=false;//LSTM Write gate. For LSTM, individual gates can be frozen separately
        frozen_fg=false;//LSTM Forget gate. For LSTM, individual gates can be frozen separately
	trait_id=1;
	override=false;
        std::cout<<"ERROR in nnode.cpp constructor. This should not have been called "<<std::endl;
        exit(0);
}

NNode::NNode(nodetype ntype,int nodeid, nodeplace placement, bool freeze, bool freeze_ip, bool freeze_rd, bool freeze_wr, bool freeze_fg, int nactivation_type) {
	active_in_flag=false;
	active_out_flag=false;
	visited=false;
	activesum=0;
	activesum_rd=0;
	activesum_wr=0;
	activesum_fg=0;
	activation=0;
        lstm_cell_state = 0.0;
	output=0;
	last_activation=0;
	last_activation2=0;
	type=ntype; //NEURON or SENSOR OR (LSTM NEURON) type
	activation_count=0; //Inactive upon creation
	node_id=nodeid;
	nodetrait=0;
	gen_node_label=placement;
        activation_type = nactivation_type;//(0-Identity, 1-Multiply, 2-Sigmoid, 3-Tanh, 4-ReLu) Evolving LSTM with Tensorflow

	if (gen_node_label == OUTPUT) {//ftype is ignored if the node is LSTM
                ftype=TANH; //Output nodes have tanh non-linearity
        }
        else {
                ftype=RELU; //Hidden nodes have relu non-linearity
        }
        dup=0;
	analogue=0;
	frozen=freeze;
        frozen_ip=freeze_ip;//LSTM Input data. For LSTM, individual gates can be frozen separately
        frozen_rd=freeze_rd;//LSTM Read gate. For LSTM, individual gates can be frozen separately
        frozen_wr=freeze_wr;//LSTM Write gate. For LSTM, individual gates can be frozen separately
        frozen_fg=freeze_fg;//LSTM Forget gate. For LSTM, individual gates can be frozen separately
	trait_id=1;
	override=false;
}

NNode::NNode(NNode *n,Trait *t) {
	active_in_flag=false;
	active_out_flag=false;
	visited=false;
	activation=0;
	activesum=0;
	activesum_rd=0;
	activesum_wr=0;
	activesum_fg=0;
        lstm_cell_state = 0.0;
	output=0;
	last_activation=0;
	last_activation2=0;
	type=n->type; //NEURON or SENSOR OR (LSTM NEURON) type
	activation_count=0; //Inactive upon creation
	node_id=n->node_id;
	nodetrait=0;
	gen_node_label=n->gen_node_label;
        activation_type = n->activation_type;//(0-Identity, 1-Multiply, 2-Sigmoid, 3-Tanh, 4-ReLu) Evolving LSTM with Tensorflow

	if (gen_node_label == OUTPUT) {//ftype is ignored if the node is LSTM
                ftype=TANH; //Output nodes have tanh non-linearity
        }
        else {
                ftype=RELU; //Hidden nodes have relu non-linearity
        }
	dup=0;
	analogue=0;
	nodetrait=t;
	frozen=n->frozen;
        frozen_ip=n->frozen_ip;//LSTM Input data. For LSTM, individual gates can be frozen separately
        frozen_rd=n->frozen_rd;//LSTM Read gate. For LSTM, individual gates can be frozen separately
        frozen_wr=n->frozen_wr;//LSTM Write gate. For LSTM, individual gates can be frozen separately
        frozen_fg=n->frozen_fg;//LSTM Forget gate. For LSTM, individual gates can be frozen separately
	if (t!=0)
		trait_id=t->trait_id;
	else trait_id=1;
	override=false;
}

NNode::NNode (const char *argline, std::vector<Trait*> &traits) {
	int traitnum;
	std::vector<Trait*>::iterator curtrait;

        activation=0;
        visited=false;
	activesum=0;
	activesum_rd=0;
	activesum_wr=0;
	activesum_fg=0;
        lstm_cell_state = 0.0;
	last_activation=0;
	last_activation2=0;
	activation_count=0; //Inactive upon creation

        std::stringstream ss(argline);
	//char curword[128];
	//char delimiters[] = " \n";
	//int curwordnum = 0;

	//Get the node parameters
	//strcpy(curword, NEAT::getUnit(argline, curwordnum++, delimiters));
	//node_id = atoi(curword);
	//strcpy(curword, NEAT::getUnit(argline, curwordnum++, delimiters));
	//traitnum = atoi(curword);
	//strcpy(curword, NEAT::getUnit(argline, curwordnum++, delimiters));
	//type = (nodetype)atoi(curword);
	//strcpy(curword, NEAT::getUnit(argline, curwordnum++, delimiters));
	//gen_node_label = (nodeplace)atoi(curword);

        int nodety, nodepl,nactivation_type;
        ss >> node_id >> nactivation_type >> nodety >> nodepl;
        type = (nodetype)nodety;
        activation_type = nactivation_type;//(0-Identity, 1-Multiply, 2-Sigmoid, 3-Tanh, 4-ReLu) Evolving LSTM with Tensorflow

        gen_node_label = (nodeplace)nodepl;

        //The following lines to assign ftype have been missing from the code previously. 
        //Could be a major bug during testing winners
	if (gen_node_label == OUTPUT) {//ftype is ignored if the node is LSTM
                ftype=TANH; //Output nodes have tanh non-linearity
        }
        else {
                ftype=RELU; //Hidden nodes have relu non-linearity
        }

	// Get the Sensor Identifier and Parameter String
	// mySensor = SensorRegistry::getSensor(id, param);
	frozen=false;  //TODO: Maybe change
        frozen_ip=false;//LSTM Input data. For LSTM, individual gates can be frozen separately
        frozen_rd=false;//LSTM Read gate. For LSTM, individual gates can be frozen separately
        frozen_wr=false;//LSTM Write gate. For LSTM, individual gates can be frozen separately
        frozen_fg=false;//LSTM Forget gate. For LSTM, individual gates can be frozen separately

	//Get a pointer to the trait this node points to
	if (traitnum==0) nodetrait=0;
	else {
		curtrait=traits.begin();
		while(((*curtrait)->trait_id)!=traitnum)
			++curtrait;
		nodetrait=(*curtrait);
		trait_id=nodetrait->trait_id;
	}

	override=false;
}

// This one might be incomplete
NNode::NNode (const NNode& nnode)
{
	active_in_flag = nnode.active_in_flag;
	active_out_flag = nnode.active_out_flag;
	visited = nnode.visited;
	activesum = nnode.activesum;
	activesum_rd = nnode.activesum_rd;
    	activesum_wr = nnode.activesum_wr;
	activesum_fg = nnode.activesum_fg;
	activation = nnode.activation;
        lstm_cell_state = nnode.lstm_cell_state;
	output = nnode.output;
	last_activation = nnode.last_activation;
	last_activation2 = nnode.last_activation2;
	type = nnode.type; //NEURON or SENSOR OR (LSTM NEURON) type
	activation_count = nnode.activation_count; //Inactive upon creation
	node_id = nnode.node_id;
	activation_type = nnode.activation_type;//(0-Identity, 1-Multiply, 2-Sigmoid, 3-Tanh, 4-ReLu) Evolving LSTM with Tensorflow

        ftype = nnode.ftype;
	nodetrait = nnode.nodetrait;
	gen_node_label = nnode.gen_node_label;
	dup = nnode.dup;
	analogue = nnode.dup;
	frozen = nnode.frozen;
        frozen_ip=nnode.frozen_ip;//LSTM Input data. For LSTM, individual gates can be frozen separately
        frozen_rd=nnode.frozen_rd;//LSTM Read gate. For LSTM, individual gates can be frozen separately
        frozen_wr=nnode.frozen_wr;//LSTM Write gate. For LSTM, individual gates can be frozen separately
        frozen_fg=nnode.frozen_fg;//LSTM Forget gate. For LSTM, individual gates can be frozen separately
	trait_id = nnode.trait_id;
	override = nnode.override;
}

NNode::~NNode() {
	std::vector<Link*>::iterator curlink;

	//Kill off all incoming links
	for(curlink=incoming.begin();curlink!=incoming.end();++curlink) {
		delete (*curlink);
	}
	//if (nodetrait!=0) delete nodetrait;
}

//Returns the type of the node, NEURON or SENSOR OR (LSTM NEURON)
const nodetype NNode::get_type() {
	return type;
}

//Allows alteration between NEURON and SENSOR.  Returns its argument
nodetype NNode::set_type(nodetype newtype) {
	type=newtype;
	return newtype;
}

//If the node is a SENSOR, returns true and loads the value
bool NNode::sensor_load(double value) {
	if (type==SENSOR) {

		//Time delay memory
		last_activation2=last_activation;
		last_activation=activation;

		activation_count++;  //Puts sensor into next time-step
		activation=value;
                active_in_flag = true; //Aditya: For the first time activate is called after network flush, this prevents output getting activated twice with the same input value
                active_out_flag = true; //Aditya: For the first time activate is called after network flush, this prevents output getting activated twice with the same input value
		return true;
	}
	else return false;
}

// Note: NEAT keeps track of which links are recurrent and which
// are not even though this is unnecessary for activation.
// It is useful to do so for 2 other reasons: 
// 1. It makes networks visualization of recurrent networks possible
// 2. It allows genetic control of the proportion of connections
//    that may become recurrent

//// Add an incoming connection a node
//void NNode::add_incoming(NNode *feednode,double weight,bool recur) {
//	Link *newlink=new Link(weight,feednode,this,recur);
//	incoming.push_back(newlink);
//	(feednode->outgoing).push_back(newlink);
//}
//
//// Nonrecurrent version
//void NNode::add_incoming(NNode *feednode,double weight) {
//	Link *newlink=new Link(weight,feednode,this,false);
//	incoming.push_back(newlink);
//	(feednode->outgoing).push_back(newlink);
//}

// Return activation currently in node, if it has been activated
double NNode::get_active_out() {
	if (activation_count>0)
		return activation;
	else return 0.0;
}

// Return activation currently in node from PREVIOUS (time-delayed) time step,
// if there is one
double NNode::get_active_out_td() {
	if (activation_count>1)
		return last_activation;
	else return 0.0;
}

// This recursively flushes everything leading into and including this NNode, including recurrencies
void NNode::flushback() {
	std::vector<Link*>::iterator curlink;
	
        visited=false;

	//A sensor should not flush black
	if (type!=SENSOR) {

		if (activation_count>0) {
			activation_count=0;
			activation=0;
			last_activation=0;
			last_activation2=0;
                        active_in_flag = false;      //ADITYA: BUG?? Currently, once active, a node remains active despite the network getting flushed. Can lead to network outputs getting activated with incoming 0. To fix this, reset the active_flag. This ensures that network-flush resets all active paths to outputs.
                        active_out_flag = false;      //ADITYA: BUG?? Currently, once active, a node remains active despite the network getting flushed. Can lead to network outputs getting activated with incoming 0. To fix this, reset the active_flag. This ensures that network-flush resets all active paths to outputs.
                        lstm_cell_state = 0.0;  //Reset LSTM cell value
		}

		//Flush back recursively
		for(curlink=incoming.begin();curlink!=incoming.end();++curlink) {
			//Flush the link itself (For future learning parameters possibility) 
			(*curlink)->added_weight=0;
			if ((((*curlink)->in_node)->activation_count>0))
				((*curlink)->in_node)->flushback();
		}
	}
	else {
		//Flush the SENSOR
		activation_count=0;
		activation=0;
		last_activation=0;
		last_activation2=0;

	}

}

// This recursively checks everything leading into and including this NNode, 
// including recurrencies
// Useful for debugging
void NNode::flushback_check(std::vector<NNode*> &seenlist) {
	std::vector<Link*>::iterator curlink;
	//int pause;
	std::vector<Link*> innodes=incoming;
	std::vector<NNode*>::iterator location;

	if (!(type==SENSOR)) {


		//std::cout<<"ALERT: "<<this<<" has activation count "<<activation_count<<std::endl;
		//std::cout<<"ALERT: "<<this<<" has activation  "<<activation<<std::endl;
		//std::cout<<"ALERT: "<<this<<" has last_activation  "<<last_activation<<std::endl;
		//std::cout<<"ALERT: "<<this<<" has last_activation2  "<<last_activation2<<std::endl;

		if (activation_count>0) {
			std::cout<<"ALERT: "<<this<<" has activation count "<<activation_count<<std::endl;
		}

		if (activation>0) {
			std::cout<<"ALERT: "<<this<<" has activation  "<<activation<<std::endl;
		}

		if (last_activation>0) {
			std::cout<<"ALERT: "<<this<<" has last_activation  "<<last_activation<<std::endl;
		}

		if (last_activation2>0) {
			std::cout<<"ALERT: "<<this<<" has last_activation2  "<<last_activation2<<std::endl;
		}

		for(curlink=innodes.begin();curlink!=innodes.end();++curlink) {
            location = std::find(seenlist.begin(),seenlist.end(),((*curlink)->in_node));
			if (location==seenlist.end()) {
				seenlist.push_back((*curlink)->in_node);
				((*curlink)->in_node)->flushback_check(seenlist);
			}
		}

	}
	else {
		//Flush_check the SENSOR


		std::cout<<"sALERT: "<<this<<" has activation count "<<activation_count<<std::endl;
		std::cout<<"sALERT: "<<this<<" has activation  "<<activation<<std::endl;
		std::cout<<"sALERT: "<<this<<" has last_activation  "<<last_activation<<std::endl;
		std::cout<<"sALERT: "<<this<<" has last_activation2  "<<last_activation2<<std::endl;


		if (activation_count>0) {
			std::cout<<"ALERT: "<<this<<" has activation count "<<activation_count<<std::endl;
		}

		if (activation>0) {
			std::cout<<"ALERT: "<<this<<" has activation  "<<activation<<std::endl;
		}

		if (last_activation>0) {
			std::cout<<"ALERT: "<<this<<" has last_activation  "<<last_activation<<std::endl;
		}

		if (last_activation2>0) {
			std::cout<<"ALERT: "<<this<<" has last_activation2  "<<last_activation2<<std::endl;
		}

	}

}

// Reserved for future system expansion
void NNode::derive_trait(Trait *curtrait) {

	if (curtrait!=0) {
		for (int count=0;count<NEAT::num_trait_params;count++)
			params[count]=(curtrait->params)[count];
	}
	else {
		for (int count=0;count<NEAT::num_trait_params;count++)
			params[count]=0;
	}

	if (curtrait!=0)
		trait_id=curtrait->trait_id;
	else trait_id=1;

}

// Returns the gene that created the node
NNode *NNode::get_analogue() {
	return analogue;
}

// Force an output value on the node
void NNode::override_output(double new_output) {
	override_value=new_output;
	override=true;
}

// Tell whether node has been overridden
bool NNode::overridden() {
	return override;
}

// Set activation to the override value and turn off override
void NNode::activate_override() {
	activation=override_value;
	override=false;
}


void NNode::print_to_file(std::ofstream &outFile) {
  outFile<<"node "<<node_id<<" ";
//if (nodetrait!=0) outFile<<nodetrait->trait_id<<" ";
//else outFile<<"0 ";
  outFile<<activation_type<<" ";//(0-Identity, 1-Multiply, 2-Sigmoid, 3-Tanh, 4-ReLu) Evolving LSTM with Tensorflow

  outFile<<type<<" ";
  outFile<<gen_node_label<<std::endl;
}


void NNode::print_to_file(std::ostream &outFile) {
	//outFile<<"node "<<node_id<<" ";
	//if (nodetrait!=0) outFile<<nodetrait->trait_id<<" ";
	//else outFile<<"0 ";
	//outFile<<type<<" ";
	//outFile<<gen_node_label<<std::endl;
        std::cout<<"Should not be HERE in node print to file"<<std::endl;
        exit(0);

	char tempbuf[128];
	sprintf(tempbuf, "node %d ", node_id);
	outFile << tempbuf;

	if (nodetrait != 0) {
		char tempbuf2[128];
		sprintf(tempbuf2, "%d ", nodetrait->trait_id);
		outFile << tempbuf2;
	}
	else outFile << "0 ";

	char tempbuf2[128];
	sprintf(tempbuf2, "%d %d\n", type, gen_node_label);
	outFile << tempbuf2;
}

//Find the greatest depth starting from this neuron at depth d
int NNode::depth(int d, Network *mynet) {
  std::vector<Link*> innodes=this->incoming;
  std::vector<Link*>::iterator curlink;
  int cur_depth; //The depth of the current node
  int max=d; //The max depth

  if (d>100) {
    //std::cout<<mynet->genotype<<std::endl;
    //std::cout<<"** DEPTH NOT DETERMINED FOR NETWORK WITH LOOP"<<std::endl;
    return 10;
  }

  //Base Case
  if ((this->type)==SENSOR)
    return d;
  //Recursion
  else {

    for(curlink=innodes.begin();curlink!=innodes.end();++curlink) {
      cur_depth=((*curlink)->in_node)->depth(d+1,mynet);
      if (cur_depth>max) max=cur_depth;
    }
  
    return max;

  } //end else

}
