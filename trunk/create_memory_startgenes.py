"""
Creates memory_startgenes file for NEAT
python create_memory_startgenes 717 100 
arg1: Number of Network Inputs (Excluding the bias. First input is always bias by default)
arg2: Number of outputs
arg3: Number of hidden units (not supported currently)
Output is a fully/partially connected network dumped in memory_startgenes file
"""
import sys
import random 

if __name__ == '__main__':


        f_wr = open("memory_startgenes", "w")
        num_inputs = int(sys.argv[1])
        num_outputs = int(sys.argv[2])
        num_genes = (num_inputs+1)*num_outputs #Fully connected network with a bias input
        lines = []

        #Start writing to the file
        lines.append('genomestart 1\n')
        lines.append('trait 1 0.0970703 0.00280173 0 0 0.00966515 0 0 0\n') 
        lines.append('trait 2 0.179688 0 0.0698419 0.00760081 0.015101 0.000960449 0.00907916 0.0111408\n') 
        lines.append('trait 3 0.43125 0.00691724 0 0.0792075 0 0 0.020273 0\n')
        
        input_node_id = []
        lines.append('node 1 0 1 3\n') #Bias node
        input_node_id.append(1) #Bias node id is always 1

        for i in xrange(0,num_inputs):
                temp_str = 'node '+str(i+2)+' 0 1 1\n'
                lines.append(temp_str) #Input node
                input_node_id.append(i+2)
        
        output_node_id = []
        for i in xrange(0,num_outputs):
                temp_str = 'node '+str(i+num_inputs+2)+' 0 0 2\n'
                lines.append(temp_str) #Output node
                output_node_id.append(i+num_inputs+2)

        gene_id = 1
        for i in xrange(0,len(output_node_id)):
                for j in xrange(0,len(input_node_id)):
                        from_node_id = str(input_node_id[j])
                        to_node_id = str(output_node_id[i])
                        if i == j: #To save memory, have only one link to the output node 
                                if (random.random() > 0.00): #Controls the sparsity of the network (Sparse networks are faster to run)
                                        weight = str(1.0)
                                        enable_flag = str(1)
                                else:
                                        weight = str(0.0)
                                        enable_flag = str(0)

                                recur_flag = str(0)
                                mutation_number = str(0.0)
                                temp_str = 'gene 3 '+from_node_id+' '+to_node_id+' '+weight+' '+recur_flag+' '+str(gene_id)+' '+mutation_number+' '+enable_flag+'\n'
                                lines.append(temp_str)
                                gene_id = gene_id + 1
        lines.append('genomeend 1\n')
        f_wr.writelines(lines)
        f_wr.close()


