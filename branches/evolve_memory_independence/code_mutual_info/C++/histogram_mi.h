#include <cstring>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <iterator>     // std::ostream_iterator
#include <iostream>
#include <algorithm>    // std::copy
#include <vector>

using namespace std;
int toBin(double v, int num_bin) {
	if(v>0.999) {
                v=0.999;
        }
        else {
                v =v+0.0001;//To avoid spurious error due to inaccuracy of double
        }
        return (int)((v)*(double)num_bin);//type casting to integer floors the number
}

double BintoValue(int v, int num_bin) {
        return (double)(((double)(v)/(double)num_bin)+0.0001);//type casting to integer floors the number
}

double log_anybase(double value, int base) {
        return (double) log(value)/log((double) base);
}


class histogram
{
        std::vector<int> bin_count;
        std::vector<int> num_bins;
        std::vector<double> granularity;
        int total;
  
        public:
        void initialize1d(std::vector<double> x, int x_num_bins) {
                num_bins.push_back(x_num_bins);
        	granularity.push_back((double)(1/num_bins[0]));
        	for (int i=0; i<num_bins[0]; i++) { 
                        bin_count.push_back(0);  //Initializing
                }
                for (int i=0; i<x.size(); i++) {
                        bin_count[toBin(x[i], num_bins[0])]++;//Incrementing bin count
                }
        	total=x.size();
        }
        void initialize2d(std::vector<double> x, std::vector<double> y,int x_num_bins, int y_num_bins) {
        	num_bins.push_back(x_num_bins);
        	num_bins.push_back(y_num_bins);
                granularity.push_back(1.0/(double)x_num_bins);
                granularity.push_back(1.0/(double)y_num_bins);
        	for (int i=0; i<num_bins[0]*num_bins[1]; i++) { 
                        bin_count.push_back(0);  //Initializing
                }
        	for(int i=0;i<x.size();i++) {
        		bin_count[toBin2d(x[i],y[i])]++;
        	}
        	total=x.size();
        }
        int toBin2d(double x,double y) {
        	if(x>0.999) x=0.999;
        	if(y>0.999) y=0.999;
        	return toBin(x, num_bins[0])*num_bins[1]+toBin(y, num_bins[1]);
        }
        
        double probability(double x,double y) {
        	return ((double)bin_count[toBin2d(x,y)])/((double)total);///(granularity[0]*granularity[1]);
        }
        
        //int toBin(double v, int num_bin) {
        //	if(v==1.0) {
        //                v=0.999;
        //        }
        //        else {
        //                v =v+0.0001;//To avoid spurious error due to inaccuracy of double
        //        }
        //        return (int)((v)*(double)num_bin);//type casting to integer floors the number
        //}
        double probability(double v) {
        	return ((double)bin_count[toBin(v, num_bins[0])])/((double)total);///(granularity[0]);
        }
        void print_bin_counts() {
                if (num_bins.size() > 1) {
                        for (int i=0; i<num_bins[0]*num_bins[1]; i++) {
                                std::cout<<bin_count[i]<<std::endl;
                        }
                        std::cout<<endl;
                }
                else {
                        for (int i=0; i<num_bins[0]; i++) {
                                std::cout<<bin_count[i]<<std::endl;
                        }
                        std::cout<<endl;
                }
        }
        
        static double standard_deviation(std::vector<double> xvals, int x_num_bins) {

        	histogram hx;// = new histogram();
        	hx.initialize1d(xvals,x_num_bins);
                //hx.print_bin_counts();
                double px; //Probability
                double mean = 0.0;
                double variance = 0.0;
                double std_dev = 0.0;
                double granularity_x = (double) (1.0/(double) x_num_bins); 

        	for(double x=0.0000001;x<1.0;x=x+granularity_x) {//Starting with x=0.0000001 to avoid double inaccuracy
                        px=hx.probability(x);
                        mean += px*x;
                }
        	for(double x=0.0000001;x<1.0;x=x+granularity_x) {//Starting with x=0.0000001 to avoid double inaccuracy
                        px=hx.probability(x);
                        variance += px*(pow((x-mean),2));
                }
                if (variance > 0.0) {
                        std_dev = sqrt(variance);
                        //std::cout<<" Standard Deviation: "<<std_dev<<" "<<sqrt(0.021)<<endl;
                }
                else {
                        std_dev = 0.0;
                }
                return std_dev;
        }
        
        static double entropy(std::vector<double> xvals, int x_num_bins) {
        	histogram hx;// = new histogram();
        	hx.initialize1d(xvals,x_num_bins);
                double granularity_x = (double) (1.0/(double) x_num_bins); 
                double entropy=0.0;
        	for(double x=0.0000001;x<1.0;x=x+granularity_x) {//Starting with x=0.0000001 to avoid double inaccuracy
                        double px=hx.probability(x);
                        if (px == 1) {//Diversity experiment (to ensure neuron value is not fixed)
                                entropy = 0.0;
                        }
                        else if(px!=0.0) {
        			entropy+=(-1.0)*px*(log_anybase(px, x_num_bins));
                        }
                }
                if (entropy > 1.0) {
                        std::cout<<" Entropy IS MORE THAN : 1.0 "<<entropy<<std::endl;//Uniform distb shd have max = 1.0
                        exit(0);
                }
        	return entropy;		

        }
        static double mutual_inf(std::vector<double> xvals, std::vector<double> yvals, int x_num_bins, int y_num_bins) {
        	histogram hxy;// = new histogram();
        	histogram hx;// = new histogram();
        	histogram hy;// = new histogram();
        	hx.initialize1d(xvals,x_num_bins);
        	hy.initialize1d(yvals,y_num_bins);
        	hxy.initialize2d(xvals,yvals,x_num_bins, y_num_bins);
                
                //hx.print_bin_counts();
                //hy.print_bin_counts();
                //hxy.print_bin_counts();
                
                double granularity_x = (double) (1.0/(double) x_num_bins); 
                double granularity_y = (double) (1.0/(double) y_num_bins); 
                
                double mi=0.0;
                double px, py, pxy, x, y;
        	for(x=0.0000001;x<1.0;x=x+granularity_x) {//Starting with x=0.0000001 to avoid double inaccuracy
                        px=hx.probability(x);
        		for(y=0.0000001;y<1.0;y=y+granularity_y) {
        			pxy=hxy.probability(x,y);
        			py=hy.probability(y);
                                if (px == 1 || py == 1) {//Diversity experiment (to ensure neuron value is not fixed)
                                        mi = 1;
                                }
                                else if(px!=0.0 && py!=0.0 && pxy !=0.0)
        				//mi+=pxy*((double)log2((double)(pxy/(px*py))));//*(0.05f*0.05f);
        				mi+=pxy*(log_anybase((double)(pxy/(px*py)), y_num_bins));//Assumption: x_num_bins = y_num_bins;
        		}
        	}
                if (mi > 1.0) {
                        std::cout<<" MUTUAL INFORMATION IS MORE THAN : 1.0 "<<mi<<std::endl;//Assuming x_num_bins = y_num_bins
                        exit(0);
                }
                //Scale mi between 0-1 by dividing by its maximum value
                //mi = mi/log2((double)y_num_bins); //Assuming x_num_bins = y_num_bins
        	return mi;		
        }

};
	
