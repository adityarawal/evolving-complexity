from pylab import *
import random, math
import histogram_probability
import numpy

def compute_mi(x, y, x_num_bins, y_num_bins):
        hx = histogram_probability.histogram_probability()
        hy = histogram_probability.histogram_probability()
        hxy = histogram_probability.histogram_probability()
        
        
        hx.initialize1d(x, x_num_bins)
        hy.initialize1d(y, y_num_bins)
        hxy.initialize2d(x,y,x_num_bins, y_num_bins);
        #Hx.print_bin_counts()
        #Hy.print_bin_counts()
        #Hxy.print_bin_counts()
        
        mi = 0.0
        
        for i in numpy.arange(0.0000001,1.0,granularity_x): #Starting with x=0.0000001 to avoid double inaccuracy
                for j in numpy.arange(0.0000001,1.0,granularity_y):
        		pxy=hxy.probability2d(i,j)
        		px=hx.probability1d(i);
        		py=hy.probability1d(j);
                        if (px == 1 or py == 1):
                                mi = 1
                                print 'px or py is one'
                                return mi
                        if(px!=0.0 and py!=0.0 and pxy !=0.0):
        			mi = mi + pxy*(math.log(1.0*(pxy/(px*py)), x_num_bins));#Assumption: x_num_bins=y_num_bins
        return mi

list_len = 10000 
x = [] #List of random numbers
y = [] #List of random numbers
z = []
x_num_bins =10
y_num_bins =10 
granularity_x =  (1.0/ x_num_bins); 
granularity_y =  (1.0/ y_num_bins); 


#for i in xrange(0,list_len):
#        x.append((random.random())) #List of random numbers
#        y.append((random.random())) #List of random numbers
#        z.append(round(random.random())) #List of random numbers
#       #x = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#        #y = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]

h1 = []
h2 = []
f1 = open("/tmp/p1.txt")
f2 = open("/tmp/p2.txt")
h1 = f1.read().strip("\n").split("\n")
h2 = f2.read().strip("\n").split("\n")
h1 = [float(i) for i in h1]
h2 = [float(i) for i in h2]
#h1 = [(float(i*j)) for i,j in zip(x,y)] 
#h2 = [(float(i*i*j*j)) for i,j in zip(x,y)] 

#from mpl_toolkits.mplot3d.axes3d import Axes3D
#
#fig = figure(figsize=(14,6))
#
## `ax` is a 3D-aware axis instance because of the projection='3d' keyword argument to add_subplot
#ax = fig.add_subplot(1, 2, 1, projection='3d')
#
#phi_m = linspace(0.001, 1, 5)
#phi_p = linspace(0.001, 1, 5)
#X,Y = meshgrid(phi_p, phi_m)
#xlabel('X',fontsize=18)
#ylabel('Y',fontsize=18)
##Z = X*Y
#Z = X/Y
#p = ax.plot_surface(X, Y, Z, rstride=4, cstride=4, linewidth=0)
##plot (range(-5,5), 5*x, label='X * Y')
##plot (range(-5,5), 5/x, label='X / Y')
#show()

#print h1
#print h2
#for i in xrange(0, list_len):
#        h1 += [(x[i] -0.8) if (x[i]-0.8 >0) else 0]
#        h2 += [(y[i] -0.8) if (y[i]-0.8 >0) else 0]
        #h1 += [float((x[i]) * (x[i]))] 
        #h2 += [float(x[i])] 
        #h2 += [float(bool(x[i]) and bool(y[i]))] 
#h1 = [ float(i) for i in x1 if i is not '']
#h2 = [ float(i) for i in x2 if i is not '']
#print h1, h2
mi = compute_mi(h1, h2, x_num_bins, y_num_bins)
print 'MI between h1 and h2: ', mi
#mi = compute_mi(h1, x, x_num_bins, y_num_bins)
#print 'MI between h1 and x: ', mi
#mi = compute_mi(h1, y, x_num_bins, y_num_bins)
#print 'MI between h1 and y: ', mi
#mi = compute_mi(h2, x, x_num_bins, y_num_bins)
#print 'MI between h2 and x: ', mi
#mi = compute_mi(h2, y, x_num_bins, y_num_bins)
#print 'MI between h2 and y: ', mi
