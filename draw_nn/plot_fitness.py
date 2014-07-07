"""
Plots the Primary and Secondary Fitness for Multiobjective NEAT and 
the primary fitness for NEAT
Currently customized for the memory recall experiment
Usage: 
python draw_neat/plot_fitness.py log/ ../evolving_complexity/log/
"""

from pylab import *
import sys, os
import fnmatch
import numpy

str1 = 'overall_average'
max_gen = 500
plot_gen = 500 
run_complete_str = ['Generation '+str(max_gen), 'WINNER IS']
max_fitness = 100 
parameter_str = 'y_x_delay'
y_x_delay = 0 
active_time_steps = 0
num_bin = 0

def mean_fitness(dirname, fitness_word_index):

        global y_x_delay 
        global active_time_steps
        global num_bin
        s1 = [max_fitness]*max_gen
        s2 = [max_fitness]*max_gen
        lists_s1 = []
        for fname in os.listdir(dirname):
            if fnmatch.fnmatch(fname, 'output*'):
                f_rd = open(dirname+'/'+fname)
                count_gen = 0
                flag_1 = 0
                flag_2 = 0
        
                for line in f_rd:
                        if parameter_str in line: #Extract parameter settings and display it in the plot  
                                words = line.strip().split()
                                y_x_delay = int(words[1])
                                active_time_steps = int(words[3])
                                num_bin = int(words[5])

                        if run_complete_str[0] in line or run_complete_str[1] in line:
                                flag_1 = 1
        
                        if str1 in line:
                                words = line.strip().split()
                                fitness1 = words[fitness_word_index] #12/5
                                s1[count_gen] = float(fitness1)
                                count_gen = count_gen + 1
                                if flag_1 == 1:
                                        flag_2 = 1
        
                        if flag_1 == 1 and flag_2 == 1:
                                flag_1 = 0 
                                flag_2 = 0
                                count_gen = 0
                                lists_s1 = lists_s1 + [s1]
                                s1 = [max_fitness]*max_gen
        
        #Find Mean of lists of varying lengths
        mean_fitness1_list = [] 
        count1 = [] 
        
        for i in xrange(max_gen):
                for j in xrange(len(lists_s1)):
                        if i < len(lists_s1[j]):
                                if i >= len(mean_fitness1_list): 
                                        mean_fitness1_list += [0]
                                        count1 += [0]
                                mean_fitness1_list[i] = mean_fitness1_list[i] + lists_s1[j][i]
                                count1[i] = count1[i] + 1  #counts the number of runs for each generation
        
        
        for i in xrange(len(mean_fitness1_list)): 
                mean_fitness1_list[i] = mean_fitness1_list[i]/count1[i]
        return mean_fitness1_list

if __name__ == "__main__":
    
        dirname = str(sys.argv[1])
        mean_fitness1_list = mean_fitness(dirname, 12) #Primary fitness in NEAT + Mutual Info
        mean_fitness2_list = mean_fitness(dirname, 15) #Secondary fitness (Mutual Info)
        t = arange(0, len(mean_fitness1_list), 1)
        plot(t[0:plot_gen], mean_fitness1_list[0:plot_gen], label='NEAT + Mutual Info Primary Fitness')
        plot(t[0:plot_gen], mean_fitness2_list[0:plot_gen], label='Secondary Fitness')
        

        if len(sys.argv) > 2: #Single Objective NEAT with/without speciation 
                dirname = str(sys.argv[2])
                mean_fitness3_list = mean_fitness(dirname, 7) #Primary fitness in NEAT + Mutual Info
                plot(t[0:plot_gen], mean_fitness3_list[0:plot_gen], label='NEAT only Fitness')


        
        legend(loc='center right', prop={'size':14})
        xlabel('Generations',fontsize=18)
        ylabel('Highest Fitness in a Generation (Mean of 30 runs)',fontsize=18)
        print y_x_delay
        title('y_x_delay = ' + str(y_x_delay) + ', active_steps = ' + str(active_time_steps) + ', num_bin = ' + str(num_bin),fontsize=16, fontweight="bold")
        tick_params(labelsize=14)
        grid(True)
        savefig("del.png")
        show()

                

