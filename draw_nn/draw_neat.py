"""
This script draws Neural Network corresponding to NEAT genome

Usage
python draw_neat.py neat_genome_file
"""


from __future__ import division
import sys
import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout

def draw_graph(neat_fname, graph, labels, graph_pos, node_id_list, node_color_list, edge_color_list, graph_layout,
               node_size=1600, node_alpha=0.3,
               node_text_size=12,
               edge_alpha=0.3, edge_tickness=1,
               edge_text_pos=0.3,
               text_font='sans-serif'):

    # create networkx graph
    G=nx.DiGraph() #Directed Graph
    #G=nx.Graph()  #Undirected Graph
    #G.add_weighted_edges_from([(1,2,0.5), (3,1,0.75)])
   
    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])

    # these are different layouts for the network you may try
    # shell seems to work best
    #if graph_layout is not None:
    #        if graph_layout == 'spring':
    #            graph_pos=nx.spring_layout(G)
    #        elif graph_layout == 'spectral':
    #            graph_pos=nx.spectral_layout(G)
    #        elif graph_layout == 'random':
    #            graph_pos=nx.random_layout(G)
    #        else:
    #            graph_pos=nx.shell_layout(G)
    graph_pos=graphviz_layout(G)
    #Delete floating (unconnected) nodes in the graph 
    i = 0
    while (i < len(node_id_list)): 
            if  (G.has_node(node_id_list[i])) is False: #If the node has no enabled genes or no connections in the graph
                    del node_id_list[i]
                    del node_color_list[i]
            else:
                    i = i + 1
    # draw graph
    nx.draw_networkx_nodes(G,graph_pos,node_id_list, node_size, 
                           node_color_list)
    nx.draw_networkx_edges(G,graph_pos, edgelist=graph, width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color_list)
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
                            font_family=text_font)

    if labels is None:
        labels = range(len(graph))

    edge_labels = dict(zip(graph, labels))
    nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels, 
                                 label_pos=edge_text_pos)

    # show graph
    plt.title(neat_fname)
    plt.show()

if __name__ == '__main__':
        node_id_list = [] 
        node_color_list =[]
        edge_color_list = []
        edge_dict = dict()
        neat_fname = str(sys.argv[1])
        f_rd = open (neat_fname, 'r')
        hidden_nodes = []
        input_nodes = []
        output_nodes = []
        bias_nodes = []
        print 'Input node is Blue'
        print 'Regular Hidden node is Gray'
        print 'LSTM Hidden node is Yellow'
        print 'Bias node is Green'
        print 'Output node is Red'
        for line in f_rd:
                if 'node' in line:#(node <node_id> <trait_number> <sensory-1/0> <hidden/input/output/bias>) 
                        words = line.strip().split()
                        node_id = str(words[1])
                        node_activation_type = str(words[2])
                        node_place = str(words[-1])
                        node_type = str(words[-2])
                        node_id_list += [node_id]
                        if node_place == '0' and node_type == '0': #Regular Hidden Node
                                hidden_nodes += [node_id]
                                if node_activation_type == '0':
                                    node_color_list += ['pink']
                                elif node_activation_type == '1':
                                    node_color_list += ['brown']
                                elif node_activation_type == '2':
                                    node_color_list += ['white']
                                elif node_activation_type == '3':
                                    node_color_list += ['gray']
                                elif node_activation_type == '4':
                                    node_color_list += ['orange']
                                else:
                                    print ('ERROR: Uncrecognized node activation type')
                                    import sys
                                    sys.exit()
                        elif node_place == '0' and node_type == '7': #LSTM Hidden Node
                                hidden_nodes += [node_id]
                                node_color_list += ['yellow']
                        elif node_place == '1': #Input node
                                input_nodes += [node_id]
                                node_color_list += ['blue']
                        elif node_place == '2': #Output node
                                output_nodes += [node_id]
                                node_color_list += ['red']
                        elif node_place == '3': #Bias node
                                bias_nodes += [node_id]
                                node_color_list += ['green']
                        else:
                                print 'ERROR: Node Type Invalid'
                                sys.exit()
                elif 'gene' in line: #(gene   <trait_number>   <from_node_id>   <to_node_id>   <connection_weight>   <recur_flag>   <gene_id>   <mutation_number>   <enable_flag>)
                        words = line.strip().split()
                        from_node = words[2]
                        to_node = words[3]
                        edge_tuple = (from_node, to_node)
                        edge_weight = (words[4])
                        gene_enable_flag = words[-1]
                        gene_recur_flag = words[5]
                        edge_color = None
                        if edge_tuple in edge_dict.keys():
                                print 'Parallel Edge', edge_tuple, 'Gene enable flag', gene_enable_flag, 'Gene recur flag', gene_recur_flag
                                continue
                        if gene_enable_flag=='1': #Store the gene only if it is enabled
                                if from_node == to_node:
                                        print 'Self Loop in Node:', from_node , to_node
                                        continue
                                if gene_recur_flag == '1':
                                        edge_color = 'red'
                                        print 'Recurrent Link in:', from_node, to_node
                                else:
                                        edge_color = 'blue'
                                edge_dict[edge_tuple] = [edge_weight, edge_color]
                                
                else:
                        continue

        """
        Determine Node positions in graph (Currently assumes only 1 hidden layer for positioning)
        """
        y_pos_offset = 2 #(Currently assumes only 1 hidden layer for positioning)
        graph_pos = dict() #This specifies the exact location of nodes
        x_max = y_max = 20
        x_min = y_min =-20
        curr_pos = [x_min, y_min]
        if (len(input_nodes) + len(bias_nodes)) > 0:
                x_pos_offset = (x_max - x_min)/(len(input_nodes) + len(bias_nodes))
                for node in input_nodes:
                        graph_pos[node] = curr_pos[:] #Copying lists by value
                        curr_pos[0] = curr_pos[0] + x_pos_offset #Nodes of input layer are laterally arranged
        
        if (len(bias_nodes)) > 0:
                curr_pos = curr_pos #Bias nodes are placed along with input nodes 
                for node in bias_nodes:
                        graph_pos[node] = curr_pos[:]
                        curr_pos[0] += x_pos_offset #Nodes of Output layer are laterally arranged
        
        if (len(hidden_nodes)) > 0:
                x_pos_offset = (x_max - x_min)/len(hidden_nodes)
                curr_pos = [x_min, y_min+y_pos_offset]
                for node in hidden_nodes:
                        graph_pos[node] = curr_pos[:]
                        curr_pos[0] += x_pos_offset #Nodes of hidden layer are laterally arranged

        if (len(output_nodes)) > 0:
                x_pos_offset = (x_max - x_min)/len(output_nodes)
                curr_pos = [x_min, y_min+2*y_pos_offset]
                for node in output_nodes:
                        graph_pos[node] = curr_pos[:]
                        curr_pos[0] += x_pos_offset #Nodes of Output layer are laterally arranged

        graph = []
        labels = []
        for edges in edge_dict:
                graph += [edges]#List tuples containing edges
                #labels += [edge_dict[edges][0]] #[List of weights] #Uncomment this to visualize weights
                edge_color_list += [edge_dict[edges][1]] #[List of weights] #Uncomment this to visualize weights


        draw_graph(neat_fname, graph, labels, graph_pos, node_id_list, node_color_list, edge_color_list, 
                                  graph_layout='shell', #None, shell, spring, spectral, random
                                  node_size=1600, 
                                  node_alpha=0.3,
                                  node_text_size=12,
                                  edge_alpha=0.3,
                                  edge_tickness=3,
                                  edge_text_pos=0.3,
                                  text_font='sans-serif')


