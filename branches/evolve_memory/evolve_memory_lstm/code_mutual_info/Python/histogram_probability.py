import random, math

def toBin(v, num_bin):
        if(v>0.999):
                v=0.999
        else:
                v =v+0.0001#To avoid spurious error due to inaccuracy of double
        return int(1.0*(v)*num_bin)#type casting to integer floors the number
        
def toBin2d(x, y, x_num_bins, y_num_bins):
        if(x>0.999):
                x=0.999
        if(y>0.999): 
                y=0.999
        return toBin(x, x_num_bins)*y_num_bins+toBin(y, y_num_bins)
                
class histogram_probability():

        def __init__(self):
                self.bin_count = []
                self.num_bins = [0, 0]
                self.total = 0; #Total number of samples

        def initialize1d(self, x, x_num_bins):
                self.num_bins[0] = x_num_bins
                self.total = len(x)
                for i in xrange(0,self.num_bins[0]): 
                        self.bin_count.append(0);  #Initializing
                
                for i in xrange(0, self.total):
                        bin_index = toBin(x[i], self.num_bins[0])
                        self.bin_count[bin_index] = self.bin_count[bin_index] + 1;#Incrementing bin count
                
                return self.bin_count        
        
        def initialize2d(self,x, y, x_num_bins, y_num_bins):
                self.num_bins[0] = x_num_bins
                self.num_bins[1] = y_num_bins
                self.total = len(x)
                for i in xrange(0,self.num_bins[0]*self.num_bins[0]): 
                        self.bin_count.append(0);  #Initializing
                
                for i in xrange(0, self.total): #Assumption: x and y have the same length
                        bin_index = toBin2d(x[i], y[i], self.num_bins[0], self.num_bins[1])
                        self.bin_count[bin_index] = self.bin_count[bin_index] + 1;#Incrementing bin count
                
                return self.bin_count        
        
        def probability1d(self,x):
                        bin_index = toBin(x, self.num_bins[0])
                	return (1.0*self.bin_count[bin_index])/(1.0*self.total)
                
        def probability2d(self,x, y):
                        bin_index = toBin2d(x, y, self.num_bins[0], self.num_bins[1])
                        return (1.0*self.bin_count[bin_index])/(1.0*self.total)#Assumption: x and y have the same length
        def print_bin_counts(self):
                        print self.bin_count
        
