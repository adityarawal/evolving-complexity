import numpy as np
import sys
import os 
min_value = 0.35
max_value = 0.75
step_size = 0.05
p2nv_path = '/scratch/cluster/aditya/memory_expt_files/p2nv_files/'
log_path  = '/scratch/cluster/aditya/memory_expt_files/log_files/'
dump_p2nv = True 
if dump_p2nv: 
    #for param1 in np.arange(0.25,max_value,max_value): 
    #    for param2 in np.arange(0.799,max_value, max_value):
    #        for param3 in np.arange(0.4,max_value,max_value):
    #            for param4 in np.arange(0.15,max_value,max_value):
    #                for param5 in np.arange(min_value,max_value,step_size): #mutate_add_node
    #                    for param6 in np.arange(min_value,max_value,step_size): #mutate_add_link
    #                        for param7 in np.arange(0.2,max_value,max_value):
    #                            for param8 in np.arange(min_value,max_value,step_size): #recur_only
    #                                for param9 in np.arange(min_value,max_value,step_size):#Mmutate_add_lstm_node
                                        for param10 in np.arange(min_value,max_value,step_size):#nw_size_cost_factor
                                            #param_str = str(param1)+'_'+str(param2)+'_'+str(param3)+'_'+str(param4)+'_'+str(param5)+'_'+str(param6)+'_'+str(param7)+'_'+str(param8)+'_'+str(param9)
                                            param_str = str(param10)
                                            p2nv_file = p2nv_path+'p2nv_'+param_str+'.ne'
                                            f_wr = open (p2nv_file, "w")
                                            lines = []
                                            print p2nv_file
                                            lines.append('trait_param_mut_prob 0.0\n')
                                            lines.append('trait_mutation_power 0.0\n')
                                            lines.append('linktrait_mut_sig 0.0\n')
                                            lines.append('nodetrait_mut_sig 0.0\n')
                                            lines.append('weight_mut_power 1.8\n')
                                            lines.append('recur_prob 0.00\n') #recur_prob is not used in the NEAT code
                                            lines.append('disjoint_coeff 1.0\n')
                                            lines.append('excess_coeff 1.0\n')
                                            lines.append('mutdiff_coeff 3.0\n')
                                            lines.append('compat_thresh 10.0\n')
                                            lines.append('age_significance 1.0\n')
                                            lines.append('survival_thresh 0.4\n')
                                            lines.append('mutate_only_prob 0.25\n')#+ str(param1)+'\n')
                                            lines.append('mutate_random_trait_prob 0.0\n')
                                            lines.append('mutate_link_trait_prob 0.0\n')
                                            lines.append('mutate_node_trait_prob 0.0\n')
                                            lines.append('mutate_link_weights_prob 0.8\n')#+ str(param2)+'\n')
                                            lines.append('mutate_toggle_enable_prob 0.4\n')#+ str(param3)+'\n')
                                            lines.append('mutate_gene_reenable_prob 0.15\n')#'+ str(param4)+'\n')
                                            lines.append('mutate_add_node_prob 0.01\n')#+ str(param5)+'\n')
                                            lines.append('mutate_add_link_prob 0.1\n')#'+ str(param6)+'\n')
                                            lines.append('interspecies_mate_rate 0.001\n')
                                            lines.append('mate_multipoint_prob 0.6\n')
                                            lines.append('mate_multipoint_avg_prob 0.4\n')
                                            lines.append('mate_singlepoint_prob 0.0\n')
                                            lines.append('mate_only_prob 0.2\n')#'+ str(param7)+'\n')
                                            lines.append('recur_only_prob 0.0\n')#'+ str(param8)+'\n')
                                            lines.append('pop_size 100 \n')
                                            lines.append('dropoff_age 150\n')
                                            lines.append('newlink_tries 100\n')
                                            lines.append('print_every 10000\n')
                                            lines.append('babies_stolen 0\n')
                                            lines.append('num_runs 1\n')
                                            lines.append('batch_size 1\n')
                                            lines.append('max_output_nodes 5\n')
                                            lines.append('frozen_startgenome 0\n')
                                            lines.append('input_sequence_len 4\n')
                                            lines.append('mutate_add_lstm_node_prob 0.00\n')#+ str(param9)+'\n')
                                            lines.append('nw_size_cost_factor '+str(param10)+'\n')
                                            
                                            f_wr.writelines(lines)
                                            f_wr.close()
                                 
                                            #Write Condor Script
                                            lines = []
                                            lines.append('Getenv = True\n')
                                            lines.append('Requirements = InMastodon\n')
                                            lines.append('Executable     = neat\n')
                                            lines.append('arguments = ' + p2nv_file +'\n') 
                                            lines.append('Log = '+log_path+'log_'+param_str+'.$(Process).txt\n')
                                            lines.append('Output  = '+log_path+'output_'+param_str+'.$(Process).txt\n')
                                            lines.append('Error  = '+log_path+'error_'+param_str+'.$(Process).txt\n')
                                            lines.append('+Group = "GRAD"\n') 
                                            lines.append('+Project = "Neural Networks"\n') 
                                            lines.append('+ProjectDescription = "Memory Evolution"\n') 
                                            lines.append('Queue 30\n')
                                            f_wr = open ('condorun', "w")
                                            f_wr.writelines(lines)
                                            f_wr.close()

                                            #Execute condor script
                                            os.system('/lusr/opt/condor/bin/condor_submit condorun')



