#include <algorithm>
#include <time.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include "../../kraskovmi.h"
#include "histogram_mi.h"

using namespace std;
vector <double> read_file(const char* fname) {
    ifstream File(fname,ios::in);//Image Data

    std::vector < double > input_data;
    std::string lineData;
    double d;

    while (getline(File, lineData)) {
            std::stringstream lineStream(lineData);
            while ((lineStream >> d)) {
                    input_data.push_back(d);
            }
    }
    File.close();
    return input_data;
}
int main() {

        vector <double> X;
        vector <double> Y;
        X = read_file("../del_X.txt"); //This data should be previously scaled to 0-1
        Y = read_file("../del_Y.txt"); //This data should be previously scaled to 0-1
        double max_X = *std::max_element(X.begin(), X.end());
        double min_X = *std::min_element(X.begin(), X.end());
        double max_Y = *std::max_element(Y.begin(), Y.end());
        double min_Y = *std::min_element(Y.begin(), Y.end());
        
        int k = 3;
        clock_t start, end;
        start = clock();
        double kraskov_mi = kraskov_mutual_information(k, X, Y);
        end = clock();
        std::cout << "Total k-NN Kraskov MI Calculation Time: "<< (double)(end-start)/CLOCKS_PER_SEC<< " seconds." << "\n";
        std::cout << " Kraskov Mutual Information: "<<kraskov_mi<<endl;
       
        int num_bins = 1000; 
        start = clock();
        double histogram_mi = histogram::mutual_inf(X, Y, num_bins, num_bins); 
        end = clock();
        std::cout << "Total Histogram MI Calculation Time: "<< (double)(end-start)/CLOCKS_PER_SEC<< " seconds." << "\n";
        std::cout << " Histogram Mutual Information: "<<histogram_mi<<endl;
        
        return 0;
}

        //std::cout<<"X max and min: "<<max_X<<" "<<min_X<<std::endl;
        //std::cout<<"Y max and min: "<<max_Y<<" "<<min_Y<<std::endl;
        //for (int i=0; i<X.size(); i++) {
        //        X[i] = (X[i]-min_X)/(max_X-min_X);
        //        Y[i] = (Y[i]-min_Y)/(max_Y-min_Y);
        //}
        //max_X = *std::max_element(X.begin(), X.end());
        //min_X = *std::min_element(X.begin(), X.end());
        //max_Y = *std::max_element(Y.begin(), Y.end());
        //min_Y = *std::min_element(Y.begin(), Y.end());
        //std::cout<<"X max and min: "<<max_X<<" "<<min_X<<std::endl;
        //std::cout<<"Y max and min: "<<max_Y<<" "<<min_Y<<std::endl;

