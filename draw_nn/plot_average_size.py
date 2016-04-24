"""
python ~/github/evolve_memory_lstm/draw_nn/plot_average_size.py /scratch/cluster/aditya/memory_expt_files/log_files
Plot average Network Size in Bar graphs
"""

from pylab import *
import sys, os
import fnmatch
import numpy


def average_network_size(fname, num_runs, search_str):
        avg_size = 0
        avg_gen = 0
        success_count = 0;
        total_runs = 0
        for run_count in np.arange(0,num_runs,1):#nw_size_cost_factor
                temp_fname = fname
                temp_fname = temp_fname+str(run_count) + '.txt'
                #Check if the file has been dumped
                if (os.path.isfile(temp_fname)): 
                        lines = os.popen("tail -20 " + temp_fname).readlines()
                        for line in lines:
                                if search_str in line:
                                        total_runs += 1
                                        success_count += 1
                                        words = line.strip().split()
                                        size = words[-1] #Last element in the line is the network size 
                                        gen = words[-7]
                                        avg_size += float(size)
                                        avg_gen += float(gen)
                                elif 'Failures' in line: 
                                        total_runs += 1
                                        words = line.strip().split()
                                        flag = int(words[1])
                                        #Success
                                        if flag == 0:
                                                total_runs -= 1 #Decrement total_runs because it will be incremented again with the winner info

        if success_count!=0:
                avg_size = avg_size/success_count
                avg_gen = avg_gen/success_count
        return (avg_size, avg_gen, success_count, total_runs) 
                

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
        min_val = 0.1
        max_val = 0.60
        step_sz = 0.05
        num_runs = 30
        search_str = "outputs 3 Size IS"
        avg_nw_size = [None]*len(np.arange(min_val,max_val,step_sz))
        avg_win_gen = [None]*len(np.arange(min_val,max_val,step_sz))
        num_success = [None]*len(np.arange(min_val,max_val,step_sz))
        total_runs = [None]*len(np.arange(min_val,max_val,step_sz))
        i = 0
        for param10 in np.arange(min_val,max_val,step_sz):#nw_size_cost_factor
            fname = dirname+'/output_'+str(param10)+'.'
            (avg_nw_size[i], avg_win_gen[i], num_success[i], total_runs[i]) = average_network_size(fname, num_runs, search_str) #Primary fitness in NEAT + Mutual Info
            avg_nw_size[i] = shorten_float(avg_nw_size[i]) #Round-off to only one decimal precision (For priniting)
            print fname, '(avg_nw_size, avg_win_gen, num_success, total_runs):',avg_nw_size[i], int(avg_win_gen[i]), num_success[i], total_runs[i]
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


        
