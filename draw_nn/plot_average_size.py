"""
Plot average Network Size in Bar graphs
"""

from pylab import *
import sys, os
import fnmatch
import numpy


def average_network_size(fname):
        temp_fname = '/tmp/del.txt'
        cmd = 'grep -i "Size IS" '+fname+' > ' + temp_fname
        os.system(cmd)
        f_rd = open(temp_fname)
        avg_size = 0
        avg_gen = 0
        count = 0;
        for line in f_rd:
                count += 1
                words = line.strip().split()
                size = words[-1] #Last element in the line is the network size 
                gen = words[-4] 
                avg_size += float(size)
                avg_gen += float(gen)
        f_rd.close()
        if count!=0:
                avg_size = avg_size/count
                avg_gen = avg_gen/count
        return (avg_size, avg_gen, count) 



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



if __name__ == "__main__":

        dirname = str(sys.argv[1])
        min_val = 0.0
        max_val = 0.2
        step_sz = 0.02
        avg_nw_size = [None]*len(np.arange(min_val,max_val,step_sz))
        avg_win_gen = [None]*len(np.arange(min_val,max_val,step_sz))
        num_finished = [None]*len(np.arange(min_val,max_val,step_sz))
        i = 0
        for param10 in np.arange(min_val,max_val,step_sz):#nw_size_cost_factor
            fname = dirname+'/output_'+str(param10)+'.*'
            (avg_nw_size[i], avg_win_gen[i], num_finished[i]) = average_network_size(fname) #Primary fitness in NEAT + Mutual Info
            print fname, '(avg_nw_size, avg_win_gen, num_finished):',avg_nw_size[i], avg_win_gen[i], num_finished[i]
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
        plot_bar_func(num_finished,min_val,max_val,step_sz, label_x, label_y, plot_title, plot_color)
        
        savefig("/u/aditya/public_html/figures/nw_cost_2.png")
        show()


        
