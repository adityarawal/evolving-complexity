import numpy as np
import sys
import os 
max_value = 0.8
step_size = 0.2
min_value = 0.1
p2nv_path = '/scratch/cluster/aditya/memory_expt_files/p2nv_files/'
log_path = '/scratch/cluster/aditya/memory_expt_files/log_files/'
dump_p2nv = True 

if dump_p2nv: 
    for param1 in np.arange(0.25,1.0,1.0): 
        for param2 in np.arange(min_value,max_value,step_size):
            for param3 in np.arange(min_value,max_value,step_size):
                for param4 in np.arange(0.15,1.0,1.0):
                    for param5 in np.arange(min_value,max_value,step_size):
                        for param6 in np.arange(min_value,max_value,step_size):
                            for param7 in np.arange(0.2,1.0,1.0):
                                for param8 in np.arange(0.0,1.0,1.0):
                                    for param9 in np.arange(min_value,max_value,step_size):
                                        param_str = str(param1)+'_'+str(param2)+'_'+str(param3)+'_'+str(param4)+'_'+str(param5)+'_'+str(param6)+'_'+str(param7)+'_'+str(param8)+'_'+str(param9)
                                        f_wr = open (p2nv_path+'p2nv_'+param_str+'.ne', "w")
                                        lines = []
                                        lines.append('trait_param_mut_prob 0.0\n')
                                        lines.append('trait_mutation_power 0.0\n')
                                        lines.append('linktrait_mut_sig 0.0\n')
                                        lines.append('nodetrait_mut_sig 0.0\n')
                                        lines.append('weight_mut_power 1.8\n')
                                        lines.append('recur_prob 0.0\n') #recur_prob is not used in the NEAT code
                                        lines.append('disjoint_coeff 1.0\n')
                                        lines.append('excess_coeff 1.0\n')
                                        lines.append('mutdiff_coeff 3.0\n')
                                        lines.append('compat_thresh 100000.0\n')
                                        lines.append('age_significance 1.0\n')
                                        lines.append('survival_thresh 0.4\n')
                                        lines.append('mutate_only_prob '+ str(param1)+'\n')
                                        lines.append('mutate_random_trait_prob 0.0\n')
                                        lines.append('mutate_link_trait_prob 0.0\n')
                                        lines.append('mutate_node_trait_prob 0.0\n')
                                        lines.append('mutate_link_weights_prob '+ str(param2)+'\n')
                                        lines.append('mutate_toggle_enable_prob '+ str(param3)+'\n')
                                        lines.append('mutate_gene_reenable_prob '+ str(param4)+'\n')
                                        lines.append('mutate_add_node_prob '+ str(param5)+'\n')
                                        lines.append('mutate_add_link_prob '+ str(param6)+'\n')
                                        lines.append('interspecies_mate_rate 0.001\n')
                                        lines.append('mate_multipoint_prob 0.6\n')
                                        lines.append('mate_multipoint_avg_prob 0.4\n')
                                        lines.append('mate_singlepoint_prob 0.0\n')
                                        lines.append('mate_only_prob '+ str(param7)+'\n')
                                        lines.append('recur_only_prob '+ str(param8)+'\n')
                                        lines.append('pop_size 100 \n')
                                        lines.append('dropoff_age 150\n')
                                        lines.append('newlink_tries 20\n')
                                        lines.append('print_every 10000\n')
                                        lines.append('babies_stolen 0\n')
                                        lines.append('num_runs 10\n')
                                        lines.append('batch_size 1\n')
                                        lines.append('max_output_nodes 1\n')
                                        lines.append('frozen_startgenome 0\n')
                                        lines.append('input_sequence_len 3\n')
                                        lines.append('mutate_add_lstm_node_prob '+ str(param9)+'\n')
                                        
                                        f_wr.writelines(lines)
                                        f_wr.close()
for param1 in np.arange(0.25,1.0,1.0): 
    for param2 in np.arange(min_value,max_value,step_size):
        for param3 in np.arange(min_value,max_value,step_size):
            for param4 in np.arange(0.15,1.0,1.0):
                for param5 in np.arange(min_value,max_value,step_size):
                    for param6 in np.arange(min_value,max_value,step_size):
                        for param7 in np.arange(0.2,1.0,1.0):
                            for param8 in np.arange(0.0,1.0,1.0):
                                for param9 in np.arange(min_value,max_value,step_size):
                                    param_str = str(param1)+'_'+str(param2)+'_'+str(param3)+'_'+str(param4)+'_'+str(param5)+'_'+str(param6)+'_'+str(param7)+'_'+str(param8)+'_'+str(param9)
                                    p2nv_file = p2nv_path+'p2nv_'+param_str+'.ne'
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
                                    lines.append('Queue 1\n')
                                    f_wr = open ('condorun', "w")
                                    f_wr.writelines(lines)
                                    f_wr.close()
                                    os.system('/lusr/opt/condor/bin/condor_submit condorun')


