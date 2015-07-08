#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <algorithm>    // std::sort
#include <vector>
#include <cmath>        // std::abs
#include <math.h> 

using namespace std;

//Use Kraskov's nearest-neighbor approach to calculate mutual information
double kraskov_mutual_information(int k, vector<double> X, vector<double> Y);

vector<double> calc_ring_blocks(int block_x_loc, int block_y_loc, int ring_num,
                                int grid_x_max, int grid_x_min, int grid_y_max, int grid_y_min);

bool twoD_nearest_neighbor(vector <double>& list_nn, vector <double>& list_nndist, vector <double> list_of_blocks, 
                                double x, double y, vector<vector<vector<double> > > twoD_grid, int k);

int total_nearby_points(double x, vector<vector<double> > oneD_grid, double dist_knn, int oneDgrid_min, int oneDgrid_max, double oneD_slicelen);

double digama ( double x );
