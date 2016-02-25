#include "kraskovmi.h"

bool myfunction (double i,double j) { return (i<j); }

bool sort_pred(const std::pair<double,int> &left, const std::pair<double,int> &right) {
        return left.first < right.first;
}

double kraskov_mutual_information(int k, const vector<double> &X, const vector<double> &Y) {
        //Create a 2D grid for the x, y points
        
        if (X.size() != Y.size()) {
            std::cout<<" Error('X and Y must contain the same number of samples')"<<std::endl;
        }
        
        int nObs = X.size();

        double max_X = 1.0; //Assumption: since sigmoid values always range between 0-1
        double min_X = 0.0; //Assumption: since sigmoid values always range between 0-1
        
        //dz = zeros(nObs, nObs);
        //dx = zeros(nObs, nObs);
        //dy = zeros(nObs, nObs);
        
        //Divide the space into a 2-D square grid of size (k/N)^(1/2)
        //Store values/pointers of X/Y values of the random variables
        double twoD_blocklen = sqrt((double)(k*10)/nObs);// Can be changed after experimentation to suit the data distribution
        int twoD_numblocks = ceil((max_X-min_X)/twoD_blocklen); 
        twoD_blocklen = (max_X-min_X)/twoD_numblocks; //Update based on the rounded off value of num blocks
        vector<vector<vector<double> > > twoD_grid(twoD_numblocks, vector<vector<double> >(twoD_numblocks));//2D grid. For each block, list of points inside the block
        std::vector<std::vector<double> > twoD_gridloc(nObs, vector<double> (2)); //For each point, returns the x, y coordinate of the block
        int twoDgrid_x_min = 0;
        int twoDgrid_x_max = twoD_numblocks-1;
        int twoDgrid_y_min = 0;
        int twoDgrid_y_max = twoD_numblocks-1;
        //vector<vector<double> > points_knn(nObs, vector<double> (2));//x, y coordinate of the k-th nearest neighbor
        vector<double> dist_knn(nObs,-1.0);

        //Create 1-D grid for x and y variables and compute number of points
        //in axis which are strictly less than dist_knn(i) apart from x(i)
        // and y(i) respectively
        double oneD_slicelen = ((double)50/nObs);// Can be changed after experimentation to suit the data distribution
        int oneD_numslice = ceil((max_X-min_X)/oneD_slicelen);
        oneD_slicelen = (max_X-min_X)/oneD_numslice; //Update based on the rounded off value of num blocks
        vector<vector<double> > oneD_x_grid(oneD_numslice);//1D grid has x locations of the points lying inside each slice of the grid
        vector<vector<double> > oneD_y_grid(oneD_numslice);//1D grid has y locations of the points lying inside each slice of the grid
        vector<double> oneD_x_gridloc(nObs);//1D x grid slice location for each point
        vector<double> oneD_y_gridloc(nObs);//1D y grid slice location for each point
        int oneDgrid_x_min = 0;
        int oneDgrid_x_max = oneD_numslice-1;
        int oneDgrid_y_min = 0;
        int oneDgrid_y_max = oneD_numslice-1;
        vector <int> nx(nObs);//Number of x points that are strictly nearer to the point of interest than the knn  
        vector <int> ny(nObs);//Number of y points that are strictly nearer to the point of interest than the knn

        clock_t start, end, start1, end1;
        for (int i=0; i<nObs; i++) {
            //if (X[i] == 1) {
            //    X[i] = X[i] - 0.0000001; //Subtract small number to make sure last block is seleced
            //}
            //if (Y[i] == 1){
            //    Y[i] = Y[i] - 0.0000001; //Subtract small number to make sure last block is selected
            //}
            //Create 2D grid
            int block_x_loc = floor((double)X[i]/twoD_blocklen); //Change ceil to floor for python/c++
            int block_y_loc = floor((double)Y[i]/twoD_blocklen); //Change ceil to floor for python/c++
            twoD_grid[block_x_loc][block_y_loc].push_back(X[i]);
            twoD_grid[block_x_loc][block_y_loc].push_back(Y[i]);
            twoD_gridloc[i][0] = block_x_loc;  
            twoD_gridloc[i][1] = block_y_loc;  
            
            //Create 1D grid
            int slice_x_loc = floor((double)X[i]/oneD_slicelen); //Change ceil to floor for python/c++
            int slice_y_loc = floor((double)Y[i]/oneD_slicelen); //Change ceil to floor for python/c++
            oneD_x_grid[slice_x_loc].push_back(X[i]);
            oneD_y_grid[slice_y_loc].push_back(Y[i]);
            oneD_x_gridloc[i] = slice_x_loc;  
            oneD_y_gridloc[i] = slice_y_loc;  
        }
        //For each point, find the k-th nearest neighbor by checking nearby grids
        for (int i=0; i<nObs; i++) {
            bool found_knn = false;
            bool found_nn = false; 
            int ring_num = 0;
            int block_x_loc = twoD_gridloc[i][0];
            int block_y_loc = twoD_gridloc[i][1]; 
            vector <double> list_nn;
            vector <double> list_nndist;
            vector <double> temp_list_nn;
            vector <double> temp_list_nndist;
            vector <double> list_of_blocks;

        
            //Find the distance to the k nearest neighbor
            while (found_knn == false) {
        
                //Step 1 - Evaluate current ring blocks and find the nearest point
                list_of_blocks = calc_ring_blocks(block_x_loc, block_y_loc, ring_num,
                                   twoDgrid_x_max, twoDgrid_x_min, twoDgrid_y_max, twoDgrid_y_min);
            
                found_nn = twoD_nearest_neighbor(temp_list_nn, temp_list_nndist, list_of_blocks, X[i], Y[i], twoD_grid,k);
                
                if(found_nn){
                    list_nn.insert(list_nn.end(), temp_list_nn.begin(), temp_list_nn.end());
                    list_nndist.insert(list_nndist.end(), temp_list_nndist.begin(), temp_list_nndist.end());
                }
                temp_list_nn.clear();
                temp_list_nndist.clear();
                
                //If nearest point found, then validate by going to Step 2.
                if ((list_nn.size()/2)>=k) {
                    
                    //Step 2 - Validate by going to the next ring.
                    ring_num = ring_num + 1;
                    list_of_blocks = calc_ring_blocks(block_x_loc, block_y_loc, ring_num,
                                   twoDgrid_x_max, twoDgrid_x_min, twoDgrid_y_max, twoDgrid_y_min);
                    
                    found_nn = twoD_nearest_neighbor(temp_list_nn, temp_list_nndist, list_of_blocks, X[i], Y[i], twoD_grid,k);
                              
                    if(found_nn) {
                        list_nn.insert(list_nn.end(), temp_list_nn.begin(), temp_list_nn.end());
                        list_nndist.insert(list_nndist.end(), temp_list_nndist.begin(), temp_list_nndist.end());
                    }
                    temp_list_nn.clear();
                    temp_list_nndist.clear();
                    
                    //Find distance to the kth nearest neighbor
                    int temp_ind;
                    double temp_dist;
                    double large_num = 100.00;
                    for (int j=0; j<k; j++) {
                            temp_ind = min_element( list_nndist.begin(), list_nndist.end() ) - list_nndist.begin();
                            temp_dist = list_nndist[temp_ind];
                            list_nndist[temp_ind] = large_num;//Replace smallest element with a large number to find the next smallest
                    }

                    int ind_kth = temp_ind;//Index of the k nearest neighbor
                    dist_knn[i] = temp_dist;

                    ////Create a vector of list_nndist indexes
                    found_knn = true;
                 
                    //********************OLD APPROACH TO FIND DISTANCE TO K-NN*********************
                    ////Create a vector of list_nndist indexes
                    //std::vector<int> ind(list_nndist.size());
                    //for (int j=0; j<ind.size(); j++) {
                    //        ind[j] = j;
                    //}

                    ////Create a pair of list_nndist and its indexes
                    //vector<pair<double,int> > nndist_ind_pair(list_nndist.size());
                    //for (int j=0; j<ind.size(); j++) {
                    //        nndist_ind_pair[j].first = list_nndist[j];
                    //        nndist_ind_pair[j].second = ind[j];
                    //}

                    //Sort in ascending order of distance
                    //std::sort(nndist_ind_pair.begin(), nndist_ind_pair.end(), sort_pred);
                    //int ind_kth = nndist_ind_pair[k-1].second;//Index of the k nearest neighbor
                    //points_knn[i][0] = list_nn[ind_kth*2];    //x location of k-th nearest point  
                    //points_knn[i][1] = list_nn[ind_kth*2+1];  //y location of k-th nearest point
                    //dist_knn[i] = list_nndist[ind_kth];
                    //********************OLD APPROACH TO FIND DISTANCE TO K-NN*********************
                }
                    
                //If no point found, then go to the next ring and go back to Step 1    
                else{     
                    ring_num = ring_num + 1;   
                }
            }
            //Find total number of points on x grid which are strictly less than k nearest distance 
            //apart from the point of interest
            nx[i] = total_nearby_points(X[i], oneD_x_grid, dist_knn[i], oneDgrid_x_min, oneDgrid_x_max, oneD_slicelen);
            
            //%Find points on y grid ....
            ny[i] = total_nearby_points(Y[i], oneD_y_grid, dist_knn[i], oneDgrid_y_min, oneDgrid_y_max, oneD_slicelen);
        }
        double mutual_info = 0.0;
        double average = 0.0;
        for (int i=0; i<nx.size(); i++) {
                average = average + digama(nx[i] + 1);
                average = average + digama(ny[i] + 1);
        }               
        average = average/nObs;
        mutual_info = digama(k) - average + digama(nObs); //psi(k) - sum(psi(nx + 1) + psi(ny + 1)) / nObs + psi(nObs)
        
        
        //If Mutual information is negative, make it zero 
        //Mutual Information can be greater than 1
        if (mutual_info < 0.0) {
                mutual_info = 0.0;
        }
        return mutual_info;
}


vector<double> calc_ring_blocks(int block_x_loc, int block_y_loc, int ring_num,
                                int grid_x_max, int grid_x_min, int grid_y_max, int grid_y_min) {
        //Returns a vector of x, y coordinate of all the blocks in a specific ring around the
        //current block
        vector <double> list_of_blocks;
        int next_block_x;
        int next_block_y;
        
        for (int i=-ring_num; i <=ring_num; i++) {
            for (int j=-ring_num; j <=ring_num; j++) {
                next_block_x = block_x_loc + i;
                next_block_y = block_y_loc + j;
                if ((i != -ring_num && i != ring_num) && (j != -ring_num && j != ring_num) 
                                || next_block_x < grid_x_min || next_block_x > grid_x_max
                                || next_block_y < grid_y_min || next_block_y > grid_y_max){
                    continue;
                }
                else{
                    list_of_blocks.push_back(next_block_x);
                    list_of_blocks.push_back(next_block_y);
                }
            }
        }
        return list_of_blocks;
}

bool twoD_nearest_neighbor(vector <double>& list_nn, vector <double>& list_nndist, vector <double> list_of_blocks, 
                           double x, double y, vector<vector<vector<double> > > twoD_grid, int k) {

        //compute distance between each sample and its k-th nearest neighbour
        //x, y - point of interest
        //list_of_blocks - find the shortest distance between point of interest and
        //the points residing inside the list_of_blocks
       
        bool found_nn; 
        int num_blocks = list_of_blocks.size()/2;//List of blocks has the x, y location of the blocks
        int block_x_loc;
        int block_y_loc;
        vector <double> points_in_block;
        double dx, dy, dz;

        //If the point of interest lies within the list_of_blocks, then exclude the point (not the
        //duplicates) while calculating the nearest neighbor
        bool flag = false;
        for (int i = 0; i<num_blocks; i++) {
            block_x_loc = list_of_blocks[i*2];
            block_y_loc = list_of_blocks[i*2+1];
            points_in_block = twoD_grid[block_x_loc][block_y_loc];
            for (int j=0; j<points_in_block.size()/2; j++) {
                dx = abs(points_in_block[j*2]-x);
                dy = abs(points_in_block[j*2+1]-y);
                if (dx == 0 && dy==0 && flag == false) {
                    flag = true;
                    continue;
                }
                else {
                        if (dx>dy) {
                                dz = dx;
                        }
                        else {
                                dz = dy;
                        }
                        list_nndist.push_back(dz);
                        list_nn.push_back(points_in_block[j*2]);  //Copy x coordinate
                        list_nn.push_back(points_in_block[j*2+1]);//Copy y coordinate
                }
            }
        }
            
        if (list_nn.empty()) { //No point found
            found_nn = false;
        }
        else { //Some points found
            found_nn = true;
        }
        return found_nn;
}

int total_nearby_points(double x, vector<vector<double> > oneD_grid, double dist_knn, int oneDgrid_min, 
                        int oneDgrid_max, double oneD_slicelen) {
        //Find points on x which are STRICTLY LESS than k nearest distance
        //apart from the point of interest
        
        int nx=0;
        int ring_num = 0;
        int init_slice_loc = floor((double)x/oneD_slicelen);
        vector <int> list_of_slices(1,init_slice_loc);
        bool found_left_side = true;
        bool found_right_side = true;
        bool flag = false;
        int slice_loc;
        vector <double> points_in_slice;
        double dx;

        //No points can be nearer than the duplicate point
        if (dist_knn == 0) {
                nx = 0;
                return nx;
        }

        while(found_left_side || found_right_side) { //Check to ensure both sides of the point are explored
            for (int i=0;i<list_of_slices.size(); i++) {
                slice_loc = list_of_slices[i];
                points_in_slice = oneD_grid[slice_loc];
                for (int j=0; j<points_in_slice.size(); j++) {
                    dx = abs(points_in_slice[j]-x);
                    if (dx == 0 && flag == false) {//Ignore the point of interest but consider other duplicate points
                        flag = true;
                        continue;
                    }
                    if (dx < dist_knn) {
                        nx = nx + 1;
                    }
                    else {
                        if (slice_loc < init_slice_loc) { //Nothing found on left side
                            found_left_side = false;
                        }
                        else if (slice_loc > init_slice_loc) { //Nothing found on right side
                            found_right_side = false;
                        }
                        else {//For points lying in the same slice as the point of interest
                            if (points_in_slice[j] < x){
                                found_left_side=false;
                            }
                            else if (points_in_slice[j] > x){
                                found_right_side = false;
                            }
                            else{
                                cout<< "Error: Duplicate point in total_nearby_points "<<x<<" "<<points_in_slice[j] <<" "<<dx<<" "<<dist_knn <<endl;
                            }
                        }
                    }
                }
            }
            
            ring_num = ring_num + 1;
            list_of_slices.clear();
            if (((init_slice_loc-ring_num) >= oneDgrid_min) && found_left_side == true) {
                list_of_slices.push_back(init_slice_loc-ring_num);
            }
            else {
                found_left_side = false;
            }
            if (((init_slice_loc+ring_num) <= oneDgrid_max) && found_right_side == true){
                list_of_slices.push_back(init_slice_loc+ring_num);
            }
            else {
                found_right_side = false;
            }
        }
        return nx;
}


double digama ( double x)

//****************************************************************************80
//
//  Purpose:
//
//    DIGAMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 June 2013
//
//  Author:
//
//    Original FORTRAN77 version by Jose Bernardo.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jose Bernardo,
//    Algorithm AS 103:
//    Psi ( Digamma ) Function,
//    Applied Statistics,
//    Volume 25, Number 3, 1976, pages 315-317.
//
//  Parameters:
//
//    Input, double X, the argument of the digamma function.
//    0 < X.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    1, X <= 0.
//
//    Output, double DIGAMA, the value of the digamma function at X.
//
{
  double euler_mascheroni = 0.57721566490153286060;
  double r;
  double value;
  double x2;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    value = 0.0;
    //*ifault = 1;
    return value;
  }
//
//  Initialize.
//
  //*ifault = 0;
  x2 = x;
  value = 0.0;
//
//  Use approximation for small argument.
//
  if ( x2 <= 0.00001 )
  {
    value = - euler_mascheroni - 1.0 / x2;
    return value;
  }
//
//  Reduce to DIGAMA(X + N).
//
  while ( x2 < 8.5 )
  {
    value = value - 1.0 / x2;
    x2 = x2 + 1.0;
  }
//
//  Use Stirling's (actually de Moivre's) expansion.
//
  r = 1.0 / x2;
  value = value + log ( x2 ) - 0.5 * r;
  r = r * r;
  value = value 
    - r * ( 1.0 / 12.0
    - r * ( 1.0 / 120.0 
    - r *   1.0 / 252.0 ) );

  return value;
}

double slow_kraskov_mutual_information(int k, const vector<double> &X, const vector<double> &Y) {

        ////List of Unique points
        //std::vector <double> X; 
        //std::vector <double> Y; 
        //
        ////Initialize by assigning the first point
        //X.push_back(W[0]);
        //Y.push_back(Z[0]);

        ////Code to generate a Unique list of 2D points from the input points (W, Z)
        ////This is necessary because Kraskov MI increases if duplicate points are present
        //for (int i=1; i <W.size(); i++) {//start comparison from the second point
        //        bool not_found = true;
        //        for (int j=0; j<X.size(); j++) {
        //                if (W[i] != X[j] || Z[i] != Y[j]) {
        //                        not_found = true;
        //                }
        //                else {
        //                        not_found = false; 
        //                        break;
        //                }
        //        }
        //        //If the point is new, then update it in the unique list 
        //        if (not_found == true) {
        //                X.push_back(W[i]);
        //                Y.push_back(Z[i]);
        //        }
        //}
        ////if ((X.size() != W.size()) || (Y.size() != Z.size())){
        ////    std::cout<<" Kraskov: Duplicate points in W, Z removed"<<std::endl;
        ////    for (int i=0; i<W.size(); i++) {
        ////            std::cout<<"W and Z: "<< W[i]<<" "<<Z[i]<<std::endl;
        ////    }
        ////    for (int i=0; i<X.size(); i++) {
        ////            std::cout<<"X and Y: "<< X[i]<<" "<<Y[i]<<std::endl;
        ////    }

        ////    exit(0);
        ////}
        
        //Old Code starts here 
        if (X.size() != Y.size()) {
            std::cout<<" Error('X and Y must contain the same number of samples')"<<std::endl;
        }
        int nObs = X.size();
        vector < vector <double> > dx (nObs, vector <double> (nObs, -1.0));
        vector < vector <double> > dy (nObs, vector <double> (nObs, -1.0));
        vector < vector <double> > dz (nObs, vector <double> (nObs, -1.0));
        vector <int> nx(nObs);//Number of x points that are strictly nearer to the point of interest than the knn  
        vector <int> ny(nObs);//Number of y points that are strictly nearer to the point of interest than the knn
        vector<double> dist_knn(nObs,-1.0);
        
        for (int i = 0; i < nObs; i++) {
            for (int j = 0; j < nObs; j++) {
                dx[i][j] = abs(X[i]- X[j]);
                dy[i][j] = abs(Y[i]- Y[j]);
                if (dx[i][j] > dy[i][j]) {
                        dz[i][j] =dx[i][j]; 
                }
                else {
                        dz[i][j] =dy[i][j]; 
                }
            }
        }
        
        for (int i = 0; i < nObs; i++) {
        
            vector <double> dxSample = dx[i];
            vector <double> dySample = dy[i];
            vector <double> dzSample = dz[i];
           
            dxSample.erase(dxSample.begin() + i);
            dySample.erase(dySample.begin() + i);
            dzSample.erase(dzSample.begin() + i); 
            
            int temp_ind;
            double temp_dist;
            double large_num = 100.00;
            for (int j=0; j<k; j++) {
                    temp_ind = min_element( dzSample.begin(), dzSample.end() ) - dzSample.begin();
                    temp_dist = dzSample[temp_ind];
                    dzSample[temp_ind] = large_num;//Replace smallest element with a large number to find the next smallest
            }

            int ind_kth = temp_ind;//Index of the k nearest neighbor
            dist_knn[i] = temp_dist;

            bool found_flag = false;
            
            //X-axis 
            std::sort (dxSample.begin(), dxSample.end(), myfunction); 
            for (int j=0; j < dxSample.size(); j++) {
                    if (dxSample[j]<dist_knn[i]) {
                            nx[i]+= 1;
                    }
                    else {
                            break;
                    }
            }
            
            //Y-axis 
            std::sort (dySample.begin(), dySample.end(), myfunction); 
            for (int j=0; j < dySample.size(); j++) {
                    if (dySample[j]<dist_knn[i]) {
                            ny[i]+= 1;
                    }
                    else {
                            break;
                    }
            }
        }
        
        double mutual_info = 0.0;
        double average = 0.0;
        for (int i=0; i<nx.size(); i++) {
                average = average + digama(nx[i] + 1);
                average = average + digama(ny[i] + 1);
        }               
        average = average/nObs;
        mutual_info = digama(k) - average + digama(nObs); //psi(k) - sum(psi(nx + 1) + psi(ny + 1)) / nObs + psi(nObs)
        
        
        //If Mutual information is negative, make it zero 
        //Mutual Information can be greater than 1
        if (mutual_info < 0.0) {
                mutual_info = 0.0;
        }
        return mutual_info;
}
