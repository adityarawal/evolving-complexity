"""
python ~/github/evolve_memory_lstm/draw_nn/plot_average_size.py /scratch/cluster/aditya/memory_expt_files/log_files
Plot average Network Size in Bar graphs
"""

from pylab import *
import sys, os
import fnmatch
import numpy
import string 


def average_network_size(fname, num_runs, num_features, task_search_str, feature_search_str):
        avg_size = 0
        avg_gen = 0
        success_count = 0;
        num_completed_runs = 0
        feature_success_list = [0]*num_features
        for run_count in np.arange(0,num_runs,1):#nw_size_cost_factor
                temp_fname = fname
                temp_fname = temp_fname+str(run_count) + '.txt'
                #Check if the file has been dumped
                if (os.path.isfile(temp_fname)): 
                        #lines = os.popen("tail -20 " + temp_fname).readlines()
                        for line in open(temp_fname, 'r'):

                                #Extract success/fail information for each feature
                                for i in xrange(0, num_features):
                                        temp_feature_str = string.replace(feature_search_str, "X", str(i+1))
                                        #temp_feature_str = list(feature_search_str)
                                        #temp_feature_str[-9] = str(i+1) #Replace X in "Num outputs X Size IS"
                                        #"".join(temp_feature_str)
                                        if (temp_feature_str in line) and ("INFOMAX WINNER" in line):
                                                feature_success_list[i] += 1
                                
                                #Extract success/fail information for the task
                                if task_search_str in line:
                                        num_completed_runs += 1
                                        success_count += 1
                                        words = line.strip().split()
                                        size = words[-1] #Last element in the line is the network size 
                                        gen = words[-7]
                                        avg_size += float(size)
                                        avg_gen += float(gen)
                                elif 'Failures' in line: 
                                        num_completed_runs += 1
                                        words = line.strip().split()
                                        flag = int(words[1])
                                        #Success
                                        if flag == 0:
                                                num_completed_runs -= 1 #Decrement num_completed_runs because it will be incremented again with the winner info

        if success_count!=0:
                avg_size = avg_size/success_count
                avg_gen = avg_gen/success_count
        return (avg_size, avg_gen, success_count, num_completed_runs, feature_success_list) 
                

def plot_bar_func (y,x_min,x_max,step_sz, label_x, label_y, plot_title, plot_color):
        N = len(y)
        x = range (N)
        ox = (np.arange(x_min,x_max,step_sz))
        bar_width = 0.5
        bar(x, y, bar_width, color=plot_color, label='Average Network size')
        title(plot_title)
        xlabel(label_x,fontsize=12)
        ylabel(label_y,fontsize=12)
        index = np.arange(N)
        xticks(index + bar_width, ox)


#Round-off to only one decimal precision (For priniting)
def shorten_float(float_num):
        str_list = [i for i in str(float_num)]
        if (len(str_list) < 3): 
                print str_list
                return float_num
        short_num = str_list[0]+str_list[1]+str_list[2]
        return (float(short_num))


if __name__ == "__main__":

        dirname = str(sys.argv[1])
        min_val = 0.45 #0.35
        max_val = 0.50 #0.75
        step_sz = 0.05
        min_fitness_thresh_value = 95.0
        max_fitness_thresh_value = 100.0
        fitness_thresh_step_size = 4.0
        num_runs = 30
        num_features = 5
        task_search_str = "TASK WINNER"
        feature_search_str = "Num outputs X Size IS"
        total_param_config =  len(np.arange(min_val,max_val,step_sz))*(len(np.arange(min_fitness_thresh_value,max_fitness_thresh_value,fitness_thresh_step_size))**4)
        avg_nw_size = [None]*total_param_config
        avg_win_gen = [None]*total_param_config
        num_success = [None]*total_param_config
        num_completed_runs = [None]*total_param_config
        i = 0
        for param10 in np.arange(min_val,max_val,step_sz):#nw_size_cost_factor
         for param11 in np.arange(min_fitness_thresh_value,max_fitness_thresh_value,fitness_thresh_step_size):
          for param12 in np.arange(min_fitness_thresh_value,max_fitness_thresh_value,fitness_thresh_step_size):
           for param13 in np.arange(min_fitness_thresh_value,max_fitness_thresh_value,fitness_thresh_step_size):
            for param14 in np.arange(min_fitness_thresh_value,max_fitness_thresh_value,fitness_thresh_step_size):
             fname = dirname+'/output_'+str(param10)+'_'+str(param11)+'_'+str(param12)+'_'+str(param13)+'_'+str(param14)+'.'
             (avg_nw_size[i], avg_win_gen[i], num_success[i], num_completed_runs[i], feature_success_list) = average_network_size(fname, num_runs, num_features, task_search_str, feature_search_str) #Primary fitness in NEAT + Mutual Info
             avg_nw_size[i] = shorten_float(avg_nw_size[i]) #Round-off to only one decimal precision (For priniting)
             print fname, '(avg_nw_size, avg_win_gen, num_success, num_completed_runs):',avg_nw_size[i], int(avg_win_gen[i]), num_success[i], num_completed_runs[i]
             for f in xrange(0, num_features): 
                     print 'feature ', str(f+1), 'success',feature_success_list[f]  
             i += 1

        sys.exit()
        figure(1)
        subplot(211)
        label_x = 'Network Size Cost Factor'
        label_y = 'Average Network Size'
        plot_title = 'Average Network Size over 30 runs v/s. Cost Factor'
        plot_color = 'blue'
        plot_bar_func(avg_nw_size,min_val,max_val,step_sz, label_x, label_y, plot_title, plot_color)
        
        subplot(212)
        label_x = 'Network Size Cost Factor'
        label_y = 'Number of times winner was found'
        plot_title = 'Number of Winning Runs (out of 30 runs) v/s. Cost Factor'
        plot_color = 'red'
        plot_bar_func(num_success,min_val,max_val,step_sz, label_x, label_y, plot_title, plot_color)
        
        savefig("/u/aditya/public_html/figures/nw_cost_2.png")
        show()


        
